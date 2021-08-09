#include"SPP.h"
#define B_MAX_SIZE 200									//矩阵长度最大值
#define L_MAX_SIZE 40
int SPP::Solve(const UserObsStore& usrStore, ostream& os)
{
	int res;
	res = Solve(usrStore);
	dumpResult(os);
	return res;
}
int SPP::Solve(const UserObsStore& usrStore)
{
	int res;
	result.Empty();
	result.epoch = usrStore.epoch;
	res = solvePosition_LS(usrStore);
	if (res == 0) {
		cerr << "SPP fault!" << endl;
		return 0;
	}
	res = modelVelocity_LS(usrStore);
	if (res == 0) {
#ifdef DEBUG
		cerr << "SPV fault!" << endl;
#endif // DEBUG
		return 0;
	}
	return 1;
}

int SPP::solvePosition_LS(const UserObsStore& usrStore)
{
	result.Empty();
	result.SYSnum[GPS] = 0;
	result.SYSnum[BDS] = 0;
	for (auto itobs:usrStore.userobsMap)
	{
		if (itobs.first.sys == GPS) result.SYSnum[GPS]++;
		else if (itobs.first.sys == BDS)result.SYSnum[BDS]++;
		else continue;
	}
	if (result.SYSnum[BDS] == 0) {													//GPS			
		return(modelPosSingleSys_LS(usrStore, GPS));
	}
	else if (result.SYSnum[GPS] == 0) {												//BDS
		return(modelPosSingleSys_LS(usrStore, BDS));
	}
	else if (result.SYSnum[BDS] != 0 && result.SYSnum[GPS] != 0) {							//双系统
		return(modelPosDoubleSys_LS(usrStore));
	}
	else {
		return 0;
	}
}

int SPP::modelPosSingleSys_LS(const UserObsStore& usrStore,SYSTEM sys)
{
	if (usrStore.size() <= 3)									//判断观测值是否达到四个
		return 0;
	result.epoch = usrStore.epoch;
	CoorTrans cTrans;														//计算工具类
	MatrixHelper m;
	TropoCorr tropoCorr;
	int no_iter = 10;														//最大迭代次数
	Triple RecX;
	Triple SatX;
	Triple RecBLH;
	Triple RecNEU;
	Triple SatNEU;
	double R0;																//几何距离
	double trop = 0;														//对流层改正数	unit:m
	double elev = 0;														//卫星高度角    unit:rad
	double recclockbias = 0;												//接收机钟差
	double V[L_MAX_SIZE];
	double PDOP;															//PDOP
	double sigma = 0;														//中误差

	for (int i = 0; i < no_iter; i++)										//开始迭代
	{
		double B[B_MAX_SIZE];												//矩阵赋值
		double L[L_MAX_SIZE] = { 0 };										//SatNum*1
		double B_T[B_MAX_SIZE];												//4*SatNum
		double Nbb[16] = { 0 }, Nbb_inv[16];								//4*4
		double B_TL[4] = { 0 };												//4*1
		double dX[4] = { 0 };												//4*1
		int SatNum = 0;														//实际可以计算的卫星数
		for (auto itobs : usrStore.userobsMap)								//B,L阵赋值
		{
			if (itobs.second.IFComb == OBSEMPTY)continue;					//剔除无组合观测值
			if (itobs.first.sys != sys) continue;							//剔除非本系统数据
			UserObsStore::UserObs& obs = itobs.second;
			SatX = obs.pvt.pos;
			if (i == 0) {
				R0 = pow(SatX[0] - RecX[0], 2) + pow(SatX[1] - RecX[1], 2)
					+ pow(SatX[2] - RecX[2], 2);							//几何距离
			}
			else {
				R0 = pow(SatX[0] - RecX[0], 2) + pow(SatX[1] - RecX[1], 2)
					+ pow(SatX[2] - RecX[2], 2);							//几何距离

				if (!cTrans.XYZToNEU(SatX, RecX, elev))						//卫星高度角计算
					continue;
				if (!cTrans.XYZToBLH(RecX, RecBLH))							//测站BLH计算
					continue;
				if (RecBLH[2] < 12000)
					trop = tropoCorr.Hopefield(elev, RecBLH[2]);			//对流层改正
				/*if (i == 6) {
					cout << itsat->first << ":  " << trop << endl;
					cout << elev * 180 / M_PI << endl;
				}*/
			}
			SatNum++;
			B[(SatNum - 1) * 4] = (RecX[0] - SatX[0]) / sqrt(R0);
			B[(SatNum - 1) * 4 + 1] = (RecX[1] - SatX[1]) / sqrt(R0);
			B[(SatNum - 1) * 4 + 2] = (RecX[2] - SatX[2]) / sqrt(R0);
			B[(SatNum - 1) * 4 + 3] = 1;
			if (sys == GPS)													//GPS L矩阵
				L[(SatNum - 1)] = obs.IFComb - sqrt(R0) + obs.pvt.clockbias*V_light
				- trop - recclockbias;
			else if (sys == BDS)
				L[(SatNum - 1)] = obs.IFComb - sqrt(R0) + obs.pvt.clockbias*V_light
				- trop - recclockbias;
			else return 0;
		}
		//最小二乘
		if (SatNum <= 3) return 0;											//实际计算的卫星数不足四个则退出
		if (!m.Matrix_Transpose(B, 4, SatNum, B_T)) {						//B_T
			return 0;
		}
		if (!m.Matrix_Mult(B_T, SatNum, 4, B, 4, SatNum, Nbb)) {			//B_T*P*B
			return 0;
		}
		if (!m.Matrix_Inv(Nbb, 4, Nbb_inv)) {								//(B_T*P*B)逆
			return 0;
		}
		if (!m.Matrix_Mult(B_T, SatNum, 4, L, 1, SatNum, B_TL)) {			//B_T*P*L
			return 0;
		}
		if (!m.Matrix_Mult(Nbb_inv, 4, 4, B_TL, 1, 4, dX)) {				//(B_T*P*B)逆*B_T*P*L
			return 0;
		}
		RecX[0] += dX[0];											//对近似值改正
		RecX[1] += dX[1];
		RecX[2] += dX[2];
		recclockbias += dX[3];

		if (sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]) < 1e-8) {	//收敛后退出
			result.SYSnum[sys] = SatNum;
			double Bx[L_MAX_SIZE] = { 0 };
			if (!m.Matrix_Mult(B, 4, SatNum, dX, 1, 4, Bx)) {
				return 0;
			}
			if (!m.Matrix_Subtract(Bx, 1, SatNum, L, 1, SatNum, V)) {		//改正数
				return 0;
			}
			for (int k = 0; k < SatNum; k++) {
				sigma += V[k] * V[k];
			}
			sigma = sqrt(sigma / (SatNum - 4.0));							//中误差
			PDOP = sqrt(Nbb_inv[0] * Nbb_inv[0]								//PDOP值
				+ Nbb_inv[5] * Nbb_inv[5] + Nbb_inv[10] * Nbb_inv[10]);
			break;
		}
		if (i == 9)
		{
			return 0;
		}
	}
	cTrans.XYZToNEU(RecX, realRecPos, RecNEU);
	result.PositionXYZ = RecX;
	result.PositionBLH[0] = RecBLH[0] / M_PI * 180;
	result.PositionBLH[1] = RecBLH[1] / M_PI * 180;
	result.PositionBLH[2] = RecBLH[2];
	result.deltaNEU = RecNEU;
	result.PositionSigma = sigma;
	result.PDOP = PDOP;

	return 1;
}

int SPP::modelPosDoubleSys_LS(const UserObsStore& usrStore)
{
	if (usrStore.size() <= 4)									//判断观测值是否达到五个
		return 0;
	result.epoch = usrStore.epoch;
	CoorTrans cTrans;														//计算工具类
	MatrixHelper m;
	TropoCorr tropoCorr;
	int no_iter = 10;														//最大迭代次数
	Triple RecX;
	Triple SatX;
	Triple RecBLH;
	Triple RecNEU;
	Triple SatNEU;

	double R0;																//几何距离
	//double traveltime = 0;												//信号传播时间	unit:s
	double trop = 0;														//对流层改正数	unit:matrix
	double elev = 0;														//卫星高度角    unit:rad
	double recclockbias_GPS = 0;											//GPS接收机钟差
	double recclockbias_BDS = 0;											//BDS接收机钟差
	double V[L_MAX_SIZE];
	double PDOP;															//PDOP
	double sigma = 0;														//中误差

	for (int i = 0; i < no_iter; i++)										//开始迭代
	{
		double B[B_MAX_SIZE];												//矩阵初始化
		double L[L_MAX_SIZE] = { 0 };										//SatNum*1
		double B_T[B_MAX_SIZE] = { 0 };										//4*SatNum
		double Nbb[25] = { 0 }, Nbb_inv[25] = { 0 };						//5*5
		double B_TL[5] = { 0 };												//5*1	
		double dX[5] = { 0 };												//5*1
		int SatNum = 0;														//实际参与计算的卫星数

		for (auto itobs : usrStore.userobsMap)								//B,L阵赋值
		{
			if (itobs.second.IFComb == OBSEMPTY)continue;							//剔除无组合观测值
			UserObsStore::UserObs& obs = itobs.second;
			SatX = obs.pvt.pos;
			if (i == 0) {
				R0 = pow(SatX[0] - RecX[0], 2) + pow(SatX[1] - RecX[1], 2)
					+ pow(SatX[2] - RecX[2], 2);							//几何距离
				//traveltime = 0.072;										//信号传播时间
			}
			else {
				R0 = pow(SatX[0] - RecX[0], 2) + pow(SatX[1] - RecX[1], 2)
					+ pow(SatX[2] - RecX[2], 2);							//几何距离

				if (!cTrans.XYZToNEU(SatX, RecX, elev))						//卫星高度角计算
					continue;
				if (!cTrans.XYZToBLH(RecX, RecBLH))							//测站BLH计算
					continue;
				if (RecBLH[2] < 12000)
					trop = tropoCorr.Hopefield(elev, RecBLH[2]);			//对流层改正
				else
					trop = 0;
			}
			SatNum++;
			B[(SatNum - 1) * 5] = (RecX[0] - SatX[0]) / sqrt(R0);
			B[(SatNum - 1) * 5 + 1] = (RecX[1] - SatX[1]) / sqrt(R0);
			B[(SatNum - 1) * 5 + 2] = (RecX[2] - SatX[2]) / sqrt(R0);
			if (itobs.first.sys == GPS) {									//GPS B,L
				B[(SatNum - 1) * 5 + 3] = 1.;
				B[(SatNum - 1) * 5 + 4] = 0.;

				L[(SatNum - 1)] = obs.IFComb - sqrt(R0) + obs.pvt.clockbias*V_light
					- trop - recclockbias_GPS;
			}
			else if (itobs.first.sys == BDS) {							//BDS B,L
				B[(SatNum - 1) * 5 + 3] = 0.;
				B[(SatNum - 1) * 5 + 4] = 1.;

				L[(SatNum - 1)] = obs.IFComb - sqrt(R0) + obs.pvt.clockbias*V_light
					- trop - recclockbias_BDS;
			}
		}
		if (SatNum <= 4) return 0;
		if (!m.Matrix_Transpose(B, 5, SatNum, B_T)) {						//B_T
			return 0;
		}
		if (!m.Matrix_Mult(B_T, SatNum, 5, B, 5, SatNum, Nbb)) {			//B_T*P*B
			return 0;
		}
		if (!m.Matrix_Inv(Nbb, 5, Nbb_inv)) {								//inv(B_T*P*B)
			return 0;
		}
		if (!m.Matrix_Mult(B_T, SatNum, 5, L, 1, SatNum, B_TL)) {			//B_T*P*L
			return 0;
		}
		if (!m.Matrix_Mult(Nbb_inv, 5, 5, B_TL, 1, 5, dX)) {				//inv(B_T*P*B)*B_T*P*L
			return 0;
		}
		RecX[0] += dX[0];											//对近似值改正
		RecX[1] += dX[1];
		RecX[2] += dX[2];
		recclockbias_GPS += dX[3];
		recclockbias_BDS += dX[4];

		if (sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]) < 1e-5) {	//收敛后退出迭代				
			double Bx[L_MAX_SIZE] = { 0 };
			if (!m.Matrix_Mult(B, 5, SatNum, dX, 1, 5, Bx)) {
				return 0;
			}
			if (!m.Matrix_Subtract(Bx, 1, SatNum, L, 1, SatNum, V)) {		//改正数
				return 0;
			}
			for (int k = 0; k < SatNum; k++) {
				sigma += V[k] * V[k];
			}
			sigma = sqrt(sigma / (SatNum - 5.0));							//中误差
			PDOP = sqrt(Nbb_inv[0] * Nbb_inv[0]								//PDOP值
				+ Nbb_inv[6] * Nbb_inv[6] + Nbb_inv[12] * Nbb_inv[12]);
			break;
		}
		if (i == 9) return 0;
	}
	cTrans.XYZToNEU(RecX, realRecPos, RecNEU);
	result.PositionXYZ = RecX;
	result.PositionBLH[0] = RecBLH[0] / M_PI * 180;
	result.PositionBLH[1] = RecBLH[1] / M_PI * 180;
	result.PositionBLH[2] = RecBLH[2];
	result.deltaNEU = RecNEU;
	result.PositionSigma = sigma;
	result.PDOP = PDOP;
	return 1;
}

int SPP::solveVel_LS(const UserObsStore& usrStore)
{
	return (modelVelocity_LS(usrStore));
}

int SPP::modelVelocity_LS(const UserObsStore& usrStore)
{
	if (usrStore.size() < 4) {												//若观测值小于4直接退出
		return 0;
	}
	if (result.PositionXYZ.norm() == 0) {								//判断是否有位置信息
		return 0;
	}
	MatrixHelper	matrix;
	int		no_iter = 1;													//最大迭代次数
	Triple	RecX = result.PositionXYZ;									//接收机位置
	Triple	SatX;
	Triple	RecV;
	Triple	SatV;
	double	clockdrift = 0;													//接收机钟速
	double	R0;																//几何距离
	double	l, m, n;														//B阵系数		
	double	V[L_MAX_SIZE];
	double	sigma = 0;

	double B[B_MAX_SIZE];												//矩阵初始化
	double L[L_MAX_SIZE] = { 0 };										//SatNum*1
	double B_T[B_MAX_SIZE] = { 0 };										//4*SatNum
	double Nbb[16] = { 0 }, Nbb_inv[16];								//4*4
	double B_TL[6] = { 0 };												//5*1	
	double dX[4] = { 0 };												//5*1
	int SatNum = 0;														//实际参与计算的卫星数
	for (auto itobs : usrStore.userobsMap)								//B L矩阵赋值
	{
		if (itobs.second.Dop[0] == OBSEMPTY)
			continue;
		SatNum++;
		UserObsStore::UserObs& obs = itobs.second;
		SatV = obs.pvt.vel;
		//cout << itsat->first << ":  " << SatV << endl;
		SatX = obs.pvt.pos;
		//cout<< itsat->first << ":  " << SatX << endl;
		R0 = pow(SatX[0] - RecX[0], 2) + pow(SatX[1] - RecX[1], 2)
			+ pow(SatX[2] - RecX[2], 2);								//几何距离
		l = (RecX[0] - SatX[0]) / sqrt(R0);
		m = (RecX[1] - SatX[1]) / sqrt(R0);
		n = (RecX[2] - SatX[2]) / sqrt(R0);
		B[(SatNum - 1) * 4] = l;										//B阵赋值
		B[(SatNum - 1) * 4 + 1] = m;
		B[(SatNum - 1) * 4 + 2] = n;
		B[(SatNum - 1) * 4 + 3] = 1;
		if (itobs.first.sys == BDS) {
			L[(SatNum - 1)] = -obs.Dop[0] * V_light / Freq_L1			//L阵赋值
				+ V_light * obs.pvt.clockdrift
				+ (l * SatV[0] + m * SatV[1] + n * SatV[2]);
		}
		else if (itobs.first.sys == GPS) {
			L[(SatNum - 1)] = -obs.Dop[0] * V_light / Freq_B1I			//L阵赋值
				+ V_light * obs.pvt.clockdrift
				+ (l * SatV[0] + m * SatV[1] + n * SatV[2]);
		}
		else
			continue;

	}
	if (SatNum < 4) return 0;											//实际计算的卫星数不足四个则退出
	if (!matrix.Matrix_Transpose(B, 4, SatNum, B_T)) {					//B_T
		return 0;
	}
	if (!matrix.Matrix_Mult(B_T, SatNum, 4, B, 4, SatNum, Nbb)) {		//B_T*P*B
		return 0;
	}
	if (!matrix.Matrix_Inv(Nbb, 4, Nbb_inv)) {							//(B_T*P*B)逆
		return 0;
	}
	if (!matrix.Matrix_Mult(B_T, SatNum, 4, L, 1, SatNum, B_TL)) {		//B_T*P*L
		return 0;
	}
	if (!matrix.Matrix_Mult(Nbb_inv, 4, 4, B_TL, 1, 4, dX)) {			//(B_T*P*B)逆*B_T*P*L
		return 0;
	}
	RecV[0] += dX[0];													//对近似值改正
	RecV[1] += dX[1];
	RecV[2] += dX[2];
	clockdrift += dX[3];

	double Bx[L_MAX_SIZE] = { 0 };
	if (!matrix.Matrix_Mult(B, 4, SatNum, dX, 1, 4, Bx)) {
		return 0;
	}
	if (!matrix.Matrix_Subtract(Bx, 1, SatNum, L, 1, SatNum, V)) {	//改正数
		return 0;
	}
	for (int k = 0; k < SatNum; k++) {
		sigma += V[k] * V[k];
	}
	sigma = sqrt(sigma / (SatNum - 4.0));							//中误差													

	result.VelocityXYZ = RecV;
	result.VelocitySigma = sigma;
	//std::cout << RecV << endl;
	return 1;
}
Triple SPP::getPos()const {
	return result.PositionXYZ;
}
void SPP::setObsID(const ObsID &id) {
	this->SPPobsID = id;
}
void SPP::dumpResult(std::ostream& os)const {
	os << std::right;									
	os.setf(ios::fixed, ios::floatfield);
	//epoch
	os << setw(4) << result.epoch.week << " " << setw(10) << setprecision(3) << result.epoch.weeksec << " ";
	//position
	os << setw(15) << setprecision(4) << result.PositionXYZ[0];		//3
	os << setw(15) << result.PositionXYZ[1];
	os << setw(15) << result.PositionXYZ[2];
	os << setw(8) << setprecision(3) << result.deltaNEU[0];		//6
	os << setw(8) << result.deltaNEU[1];
	os << setw(8) << result.deltaNEU[2];
	os << setw(14) << setprecision(8) << result.PositionBLH[0];		//9
	os << setw(14) << result.PositionBLH[1];
	os << setw(14) << result.PositionBLH[2];
	//velocity
	os << setw(10) << setprecision(3) << result.VelocityXYZ[0];		//12
	os << setw(10) << result.VelocityXYZ[1];
	os << setw(10) << result.VelocityXYZ[2];
	//精度指标
	os << setw(7) << setprecision(3) << result.PDOP;					//15
	os << setw(7) << result.PositionSigma;
	os << setw(7) << result.VelocitySigma;
	//卫星数
	os << setw(3) << setprecision(0) << result.SYSnum[GPS];				//18
	os << setw(3) << result.SYSnum[BDS];
	os << setw(3) << result.SYSnum[GPS] + result.SYSnum[BDS];
	os << std::endl;
}
void SPP::setRealPos(const Triple& pos) {
	realRecPos = pos;
}
