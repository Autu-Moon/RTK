#include "RTK.h"
const Triple RTK::BaseRec = Triple(-2267804.5263, 5009342.3723, 3220991.8632);
double RTK::RatioThreshold = 1;
int RTK::Solve(UserObsStore& base, UserObsStore& rover, bool fixed, const Triple& RoverX, const Triple& BaseX) {
	result.Empty();
	if (base.epoch != rover.epoch) {
		std::cerr << "Time doesn't match in RTK" << std::endl;
		return 2;
	}
	result.epoch = base.epoch;
	selectRoverObs(base, rover);
	//计算卫星数
	for (auto r : rover) {
		result.SatNum[r.first.sys]++;
	}
	this->doudiffMap.clear();
	Triple gpsSize = setDouDiffMap(base, rover, GPS);
	Triple bdsSize = setDouDiffMap(base, rover, BDS);
	//建模 解算
	if (!modelPosition_LS(base, rover, gpsSize + bdsSize, fixed, RoverX, BaseX)) {
		std::cerr << "RTK error in " << base.epoch << std::endl;
		return 0;
	}
#ifdef DEBUG_RTK
	dumpDouDiff(cout);
#endif // DEBUG_RTK
	doudiffMap.clear();
	return 1;
}

int	RTK::Solve(UserObsStore& base, UserObsStore& rover, ostream& os, bool fixed, const Triple& RoverX, const Triple& BaseX) {
	if (!Solve(base, rover, fixed,RoverX, BaseX)) {
		return 0;
	}
	dumpResult(os);
	return 1;
}

Triple RTK::setDouDiffMap(UserObsStore& base, UserObsStore& rover, SYSTEM sys) {
	Triple Size(0, 0, 0);
	UserObsStore::UserObs* pivotBase(0);					//Pivot - Base
	UserObsStore::UserObs* pivotRover(0);					//Pivot - Rover
	UserObsStore::UserObsMap::const_iterator normRover		//Norm - Rover
		= rover.userobsMap.cbegin();

	// if no pivot sat
	if (PivotSat[sys].isEmpty()) {
		//if no pivot sat skip this system
		if (!selectPivotSat(base, rover, sys, PivotSat[sys]))
			return Size;
	}
	pivotBase = &base.userobsMap[PivotSat[sys]];
	pivotRover = &rover.userobsMap[PivotSat[sys]];
	// if pivot-base lack of observation
	if (pivotBase->isLack() || pivotRover->isLack()) {
		//if no pivot sat skip this system
		if (!selectPivotSat(base, rover, sys, PivotSat[sys]))
			return Size;
	}
	pivotBase = &base.userobsMap[PivotSat[sys]];
	pivotRover = &rover.userobsMap[PivotSat[sys]];
	//loop by Norm satellite
	for (auto normBase : base.userobsMap) {					//Norm - Base				
		//skip pivot satellite and satellites of other system
		if (normBase.first == PivotSat[sys] || normBase.first.sys != sys) {
			continue;
		}
		normRover = rover.userobsMap.find(normBase.first);
		if (normRover == rover.userobsMap.cend())
			continue;
		//iterate by frequency 
		for (size_t i = 0; i < OBSOCCUPIED; i++) {
			//Psr
			if (!ObsData::isEmpty(normBase.second.Psr[i]) && !ObsData::isEmpty(normRover->second.Psr[i])) {
				doudiffMap[_C][sys][getFreqBand(i, sys)][normBase.first] = computeDouDiff(pivotBase->Psr[i],
					pivotRover->Psr[i], normBase.second.Psr[i], normRover->second.Psr[i]);
				Size[0]++;
			}
			//Adr
			if (!ObsData::isEmpty(normBase.second.Adr[i]) && !ObsData::isEmpty(normRover->second.Adr[i])) {
				doudiffMap[_L][sys][getFreqBand(i, sys)][normBase.first] = computeDouDiff(pivotBase->Adr[i],
					pivotRover->Adr[i], normBase.second.Adr[i], normRover->second.Adr[i]);
				Size[1]++;
			}
		}
		//store counterpart Sat PVT
	}
	return Size;
}

mym::Matrixd RTK::modelCovDx(const UserObsStore& base, const UserObsStore& rover,
	const ObsCode& code, const SYSTEM& sys, const FreqBand& freq) {
	if (doudiffMap[code][sys][freq].size() == 0) {
		return MyMatrix::Matrixd(0, 0);
	}
	std::map<SatID, double>& ddf = doudiffMap[code][sys][freq];
	size_t ddfSize = ddf.size();
	/*Cov_ND:
	*	rover_2
	*		base_2
	*			rover_3
	*				base_3
	*					......
	*						rover_pivot
	*							base_pivot
	*/
	MyMatrix::Matrixd Cov_ND(2 * (ddfSize + 1), 2 * (ddfSize + 1));
	/*C_SD:
	*	1 -1  0  0  0  0 ...
	*	0  0  1 -1  0  0 
	*	0  0  0  0  1 -1
	*
	*/
	MyMatrix::Matrixd C_SD(ddfSize + 1, 2 * (ddfSize + 1));
	/*Cov_SD:
	*	r-b_2
	*		r-b_3
	*			......
	*				r-b_pivot
	*/
	MyMatrix::Matrixd Cov_SD(ddfSize + 1, ddfSize + 1);
	/*C_DD:
	*	1  0  0  ...   0  -1
	*	0  1  0  ...   0  -1
	*	0  0  1  ...   0  -1
	*	......
	*	0  0  0  ...   1  -1
	*/
	MyMatrix::Matrixd C_DD(ddfSize, ddfSize + 1);
	/*Cov_SD:
	*	r-b_2-pivot
	*		r-b_3-pivot
	*			......
	*				r-b_n-pivot
	*/
	MyMatrix::Matrixd Cov_DD(ddfSize, ddfSize);
	Cov_ND.SetZeros();
	C_SD.SetZeros();
	Cov_SD.SetZeros();
	C_DD.SetZeros();
	Cov_DD.SetZeros();
	size_t Count = 0;
	//loop by  norm satellite
	for (auto sat : ddf) {
		Cov_ND[2 * Count][2 * Count] = weightElev(rover[sat.first].elev, code);
		Cov_ND[2 * Count + 1][2 * Count + 1] = weightElev(base[sat.first].elev, code);
		C_SD[Count][2 * Count] = 1;
		C_SD[Count][2 * Count + 1] = -1;
		C_DD[Count][Count] = 1;
		C_DD[Count][ddfSize] = -1;
		Count++;
	}
	//pivot Satellite
	Cov_ND[2 * Count][2 * Count] = weightElev(rover[PivotSat[sys]].elev, code);
	Cov_ND[2 * Count + 1][2 * Count + 1] = weightElev(base[PivotSat[sys]].elev, code);
	C_SD[Count][2 * Count] = 1;
	C_SD[Count][2 * Count + 1] = -1;
	//compute doudiff D
	Cov_SD = C_SD * Cov_ND * C_SD.Transpose();
	Cov_DD = C_DD * Cov_SD * C_DD.Transpose();
	return Cov_DD;
}

int RTK::modelPosition_LS(const UserObsStore& base, const UserObsStore& rover, const Triple& ObsCount, bool fixed,
	const Triple& RoverX, const Triple& BaseX) {
	//main rank
	/*
	*GPS L1 Psr
	*	GPS L2 Psr
	*		BDS B1I Psr
	*			BDS B3I Psr
	*				GPS L1 Adr
	*					GPS L2 Adr
	*						BDS B1I Adr
	*							BDS B3I Adr
	*/
	//x
	/*x=[x y z N_L1 N_L2 N_B1I N_B3I]*/
	//B
	/*
	*  l   m   n   0   0   0   0
	*  l   m   n   0   0   0   0
	*  l   m   n   0   0   0   0
	*  l   m   n   0   0   0   0
	*  l   m   n  l_l1 0   0   0
	*  l   m   n   0  l_l2 0   0
	*  l   m   n   0   0  l_b1 0
	*  l   m   n   0   0   0  l_b2
	*/
	size_t _r = ObsCount[0] - 3;
	size_t _t = ObsCount[1] + 3;
	size_t no_iter = 10;						//迭代次数
	MyMatrix::Matrixd B(_r + _t, _t);			//B矩阵
	MyMatrix::Matrixd L(_r + _t, 1);			//L矩阵
	MyMatrix::Matrixd P(_r + _t, _r + _t);		//观测值权阵
	MyMatrix::Matrixd Q;						//解的协方差阵
	MyMatrix::Matrixd Qbb;
	MyMatrix::Matrixd V;						//改正数
	MyMatrix::Matrixd dx;						
	MyMatrix::Matrixd N;						//模糊度固定解

	Triple RCoor = RoverX;				//流动站坐标（待估参数）
	Triple BCoor = BaseX;				//基准站坐标
	Triple pivotBase[SYScount];			//基准站算出的参考星坐标
	Triple pivotRover[SYScount];		//流动站算出的参考星坐标
	Triple normBase;					//基准站算出的非参考星坐标
	Triple normRover;					//流动站算出的非参考星坐标
	CoorTrans ctrans;
	LeastSquare LS;
	
	for (size_t sys = 0; sys < SYScount; sys++) {
		auto itbase = base.userobsMap.find(PivotSat[sys]);
		auto itrover = rover.userobsMap.find(PivotSat[sys]);
		if (itbase == base.userobsMap.end() || itrover == rover.userobsMap.end()) {
			return 0;
		}
		pivotBase[sys] = itbase->second.pvt.pos;
		pivotRover[sys] = itrover->second.pvt.pos;
	}
	
	double _nr = 0;				//卫星-测站距离
	double _pr = 0;
	double _nb = 0;
	double _pb = 0;
	double waveLength = 0;		//信号波长
	Triple Nlmn;				//B矩阵系数
	Triple Plmn;

	//model Cov Dx
	//loop by obsCode
	size_t rowCount = 0;
	size_t colCount = 0;
	for (auto obsCode : doudiffMap) {
		//loop by SYSTEM
		for (auto system : obsCode.second) {
			//loop by FreqBand
			for (auto freqband : system.second) {
				MyMatrix::Matrixd Pblock = modelCovDx(base, rover, obsCode.first, system.first, freqband.first).Inverse();
				P.SetBlock(rowCount, colCount, Pblock.Row(), Pblock.Col(), Pblock);
#ifdef DEBUG_RTK
				//std::cout << obscodeString[obsCode.first] << "  " << systemString[system.first] << "  " << freqbandString[freqband.first] << std::endl;
				//std::cout << P << std::endl;
#endif // DEBUG_RTK
				rowCount += Pblock.Row();
				colCount += Pblock.Col();
			}
		}
	}

	//Least Square Float Solution 
	for (size_t iter = 0; iter <= no_iter; iter++) {
		//initialize
		double R0 = 0;
		size_t PsrCount = 0;
		size_t AdrCount = 0;
		B.SetZeros(); 
		L.SetZeros();
		//loop by obsCode
		for(auto obsCode:doudiffMap) {
			//loop by SYSTEM
			for (auto system:obsCode.second) {
				//loop by FreqBand
				for (auto freqband: system.second) {
					//loop by satllite
					for (auto sat : freqband.second) {
						normRover = rover[sat.first].pvt.pos;
						normBase = base[sat.first].pvt.pos;
						Nlmn = computeCoeffB(RCoor, normRover);
						Plmn = computeCoeffB(RCoor, pivotRover[system.first]);
						_nr = (normRover - RCoor).norm();
						_pr = (pivotRover[system.first] - RCoor).norm();
						_nb = (normBase - BCoor).norm();
						_pb = (pivotBase[system.first] - BCoor).norm();
						//伪距观测值
						if (obsCode.first == _C) {
							B[PsrCount][0] = Nlmn[0] - Plmn[0];
							B[PsrCount][1] = Nlmn[1] - Plmn[1];
							B[PsrCount][2] = Nlmn[2] - Plmn[2];
							L[PsrCount][0] = (sat.second - ((_nr - _pr) - (_nb - _pb)));
							PsrCount++;
						}
						//载波相位观测值
						else if (obsCode.first == _L) {
							waveLength = V_light / getFreq(freqband.first);
							B[PsrCount + AdrCount][0] = Nlmn[0] - Plmn[0];
							B[PsrCount + AdrCount][1] = Nlmn[1] - Plmn[1];
							B[PsrCount + AdrCount][2] = Nlmn[2] - Plmn[2];
							B[PsrCount + AdrCount][3 + AdrCount] = waveLength;
							L[PsrCount + AdrCount][0] = (sat.second * waveLength - ((_nr - _pr) - (_nb - _pb)));
							AdrCount++;
						}
						else
							std::cerr << "RTK B fault" << std::endl;
					}
				}
			}
		}
		if (!LS.Solve(B, P, L, dx)) {
			return 0;
		}
		//更新
		RCoor[0] += dx[0][0];
		RCoor[1] += dx[1][0];
		RCoor[2] += dx[2][0];
		//收敛
		if (sqrt(dx[0][0] * dx[0][0] + dx[1][0] * dx[1][0] + dx[2][0] * dx[2][0]) < 10e-8) {
			V = B * dx - L;
			Q = (B.Transpose()*P*B).Inverse();
			Qbb = Q.Cut(0, 0, 3, 3);
			break;
		}
		if (iter == no_iter) {
			return 0;
		}
	}
	//Fixed Solution
	if (fixed) {
		result.fixed = fixAmbiguity(Q, dx, 2, N);
		if (result.fixed) {
			MyMatrix::Matrixd Xfixed(3, 1);
			MyMatrix::Matrixd Xfloat(3, 1);
			Xfloat[0][0] = RCoor[0];
			Xfloat[1][0] = RCoor[1];
			Xfloat[2][0] = RCoor[2];
			Xfixed = Xfloat - Q.Cut(0, 3, 3, Q.Col() - 3)*Q.Cut(3, 3, Q.Row() - 3, Q.Col() - 3).Inverse()*
				(dx.Cut(3, 0, dx.Row() - 3, 1) - N);
			RCoor[0] = Xfixed[0][0];
			RCoor[1] = Xfixed[1][0];
			RCoor[2] = Xfixed[2][0];
			Qbb = Qbb - Q.Cut(0, 3, 3, Q.Col() - 3)*Q.Cut(3, 3, Q.Row() - 3, Q.Col() - 3).Inverse()*
				Q.Cut(3, 0, Q.Row() - 3, 3);
		}
	}
	//结果整理
	result.PositionXYZ = RCoor;
	ctrans.XYZToNEU(result.PositionXYZ, realRecPos, result.deltaNEU);
	ctrans.XYZToNEU(result.PositionXYZ, BaseX, result.Baseline);
	result.PositionSigma = (V.Transpose() * P * V)[0][0] / _r;
	result.PDOP = Norm2({ Qbb[0][0],Qbb[1][1],Qbb[2][2] });
#ifdef DEBUG_RTK
	if (result.deltaNEU.norm() > 4) {
		std::cout << "stop" << std::endl;
	}
	/*std::ofstream fout("RTK_D.txt", ios::app);
	fout << std::right;
	fout.setf(ios::fixed);
	fout << base.epoch << "  " << ObsCount[0] << "  " << ObsCount[1] << endl;
	fout << P << endl;*/
#endif // DEBUG_RTK
	return 1;

}

bool RTK::fixAmbiguity(const mym::Matrixd& Q, const mym::Matrixd dx, int fixedNum,mym::Matrixd& N,double RatThre) {
	mym::Matrixd Qcut = Q.Cut(3, 3, Q.Row() - 3, Q.Col() - 3);
	mym::Matrixd dxcut = dx.Cut(3, 0, dx.Row() - 3, 1);
	double *ArrayQ = new double[Qcut.Col()*Qcut.Row()];
	double *ArrayFloatN = new double[dxcut.Row()];
	double *ArrayFixedN = new double[dxcut.Row()*fixedNum];
	double *ArrayFixedSigma = new double[fixedNum];
	Qcut.ExportArray(Qcut.Col()*Qcut.Row(), ArrayQ);
	dxcut.ExportArray(dxcut.Row(), ArrayFloatN);
	//LAMBDA
	lambda(dxcut.Row(), fixedNum, ArrayFloatN, ArrayQ, ArrayFixedN, ArrayFixedSigma);
	double ratio = ArrayFixedSigma[1] / ArrayFixedSigma[0];
	if (ratio < RatioThreshold) {
		return false;
	}
	N = MyMatrix::Matrixd::Zeros(dxcut.Row(), 1);
	N.ImportArray(ArrayFixedN, dxcut.Row(), 1);

	delete ArrayFloatN;
	delete ArrayQ;
	delete ArrayFixedN;
	delete ArrayFixedSigma;
	return true;
}

void RTK::selectBaseObs(UserObsStore& base, UserObsStore& rover) {
	UserObsStore::UserObsMap::const_iterator irover = rover.cbegin();
	for (UserObsStore::UserObsMap::iterator ibase=base.begin(); 
		ibase != base.end();) {
		irover = rover.userobsMap.find(ibase->first);
		if (irover == rover.cend()) {
			base.userobsMap.erase(ibase++);
		}
		else
			ibase++;
	}
}

void RTK::selectRoverObs(UserObsStore& base, UserObsStore& rover) {
	UserObsStore::UserObsMap::const_iterator ibase = base.cbegin();
	for (UserObsStore::UserObsMap::iterator irover = rover.begin();
		irover != rover.end();) {
		ibase = base.userobsMap.find(irover->first);
		if (ibase == base.cend()) {
			rover.userobsMap.erase(irover++);
		}
		else
			irover++;
	}
}

bool RTK::selectPivotSat(UserObsStore& base, UserObsStore& rover,SYSTEM sys,SatID& satid) {
	SatID sat = rover.getBest(sys);
	UserObsStore::UserObsMap::const_iterator ibase = rover.userobsMap.find(sat);
	size_t loop = 0;
	//until find the best sat both has and make sure no lack of obs
	while (ibase == rover.cend()||ibase->second.isLack()) {
		rover.userobsMap[sat].best = 0;
		sat = rover.getBest(sys);
		ibase = rover.userobsMap.find(sat);
		loop++;
		// if no best sat
		if (loop >= rover.userobsMap.size()) {
			return false;
		}
	}
	satid = sat;
	return true;
}

double RTK::computeDouDiff(double pb, double pr, double nb, double nr) {
	//return ((nr - nb) - (pr - pb));
	return ((nr - pr) - (nb - pb));
}

void RTK::dumpResult(std::ostream& os) {
	os << std::right;
	os.setf(ios::fixed, ios::floatfield);
	//epoch
	os << setw(4) << result.epoch.week << " " << setw(10) << setprecision(3) << result.epoch.weeksec << " ";
	//position
	os << setw(15) << setprecision(3) << result.PositionXYZ[0];	//3
	os << setw(15) << result.PositionXYZ[1];
	os << setw(15) << result.PositionXYZ[2];
	os << setw(9) << setprecision(4) << result.deltaNEU[0];		//6
	os << setw(9) << result.deltaNEU[1];
	os << setw(9) << result.deltaNEU[2];
	os << setw(9) << setprecision(4) << result.Baseline[0];		//6
	os << setw(9) << result.Baseline[1];
	os << setw(9) << result.Baseline[2];
	//fixed
	os << setw(3) << result.fixed;
	//os << setw(15) << setprecision(8) << result.PositionBLH[0];		//9
	//os << setw(15) << result.PositionBLH[1];
	//os << setw(15) << result.PositionBLH[2];
	//精度指标
	os << setw(10) << setprecision(7) << result.PDOP;					//15
	//os << setw(14) << result.PositionSigma;
	//卫星数
	os << setw(3) << setprecision(0) << result.SatNum[GPS];				//18
	os << setw(3) << result.SatNum[BDS];
	os << setw(3) << result.SatNum[GPS] + result.SatNum[BDS];
	os << std::endl;
}

void RTK::setObsID(const ObsID& id) {
	this->RTKobsID = id;
}

void RTK::dumpDouDiff(std::ostream& os) {
	os << result.epoch << std::endl;
	os << std::right;
	os.setf(ios::fixed, ios::floatfield);
	for (auto obsCode : doudiffMap) {
		if (obsCode.first == _C)
			os << "Psr  ";
		else
			os << "Adr  ";
		//loop by SYSTEM
		for (auto system : obsCode.second) {
			os << systemString[system.first] << "  ";
			//loop by FreqBand
			for (auto freqband : system.second) {
				os << freqbandString[freqband.first] << std::endl;
				//loop by satllite
				for (auto sat : freqband.second) {
					os << sat.first << "  " <<setprecision(9)<< std::setw(15) << sat.second << std::endl;
				}
			}
		}
	}
}

Triple RTK::computeCoeffB(const Triple& rec, const Triple& sat) {
	Triple lmn;
	double distance = (sat - rec).norm();
	lmn[0] = - (sat[0] - rec[0]) / distance;
	lmn[1] = - (sat[1] - rec[1]) / distance;
	lmn[2] = - (sat[2] - rec[2]) / distance;
	return lmn;
}

double RTK::weightElev(double elev, const ObsCode& code) {
	double alpha = 1.5;	
	double sigma_0 = 0.003*0.003;
	int PAR = 5000;
	double sigma_2 = 0;
	if (code == _C) {
		sigma_2 = sigma_0 * (1 + alpha * cos(elev) * cos(elev)) * PAR; //unit: m
	}
	else if (code == _L) {
		sigma_2 = sigma_0 * (1 + alpha * cos(elev) * cos(elev)); //unit: m
	}
	else {
		std::cerr << "ERROR in ObsCode, Only for _C&_L" << std::endl;
		return 10000000;
	}
	return sigma_2;
}

void RTK::setRealRecPos(const Triple& real) {
	this->realRecPos = real;
}