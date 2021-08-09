#include"Ephemeris.h"
int GPSEphemeris::getSatXvt(const GPSweek& t, PVT& sat)
{
	//判断星历是否过期（2h）
	if (!isValid(t))
		return 0;
	//卫星位置计算
	double A = rootA * rootA;
	//运行时间
	double tk = (t - toe).toSec();
	//toe时刻卫星平均角速度
	double n0 = sqrt(GM_GPS / (A*A*A));
	//观测时刻卫星平均角速度
	double n = n0 + deltaN;
	//平近点角
	double M = M0 + n * tk;
	M = fmod(M + 2 * M_PI, 2 * M_PI);
	//偏近点角
	double E = M;
	if (!computeEk(M, E)) {
		return 0;
	}
	//真近点角
	double f = atan2(sqrt(1 - ecc * ecc)*sin(E), cos(E) - ecc);
	//升交角距 卫星向径 轨道倾角
	double phi = f + omega;
	phi = fmod(phi, 2 * M_PI);

	double u = phi + cuc * cos(2 * phi) + cus * sin(2 * phi);
	double r = A * (1 - ecc * cos(E)) + crc * cos(2 * phi) + crs * sin(2 * phi);
	double i = i0 + i_dot * tk + cic * cos(2 * phi) + cis * sin(2 * phi);
	//升交点经度
	double Omega_k = Omega0 + (Omega_dot - Omega_e_GPS)*tk - Omega_e_GPS * toe.weeksec;
	Omega_k = fmod(Omega_k + 2 * M_PI, 2 * M_PI);
	//卫星轨道平面坐标
	double pos = r * cos(u);
	double y = r * sin(u);
	//瞬时地球坐标系坐标
	sat.pos.thearray[0] = pos * cos(Omega_k) - y * cos(i)*sin(Omega_k);
	sat.pos.thearray[1] = pos * sin(Omega_k) + y * cos(i)*cos(Omega_k);
	sat.pos.thearray[2] = y * sin(i);
	//卫星速度计算
	double E_dot = n / (1 - ecc * cos(E));
	double phi_dot = sqrt((1 + ecc) / (1 - ecc))*
		pow((cos(f / 2) / cos(E / 2)), 2)*E_dot;
	//对位置求导
	double u_dot = 2 * (cus*cos(2 * phi) - cuc * sin(2 * phi))*phi_dot + phi_dot;
	double r_dot = A * ecc*sin(E)*E_dot + 2 * (crs*cos(2 * phi) - crc * sin(2 * phi))*phi_dot;
	double I_dot = i_dot + 2 * (cis*cos(2 * phi) - cic * sin(2 * phi))*phi_dot;
	double Omega_k_dot = Omega_dot - Omega_e_GPS;
	double x_dot = r_dot * cos(u) - r * u_dot*sin(u);
	double y_dot = r_dot * sin(u) + r * u_dot*cos(u);
	//求导结果以矩阵形式展示
	double R[12] =
	{ cos(Omega_k),-sin(Omega_k)*cos(i),-(pos*sin(Omega_k) + y * cos(Omega_k)*cos(i)),y*sin(Omega_k)*sin(i),
		sin(Omega_k),cos(Omega_k)*cos(i),(pos*cos(Omega_k) - y * sin(Omega_k)*cos(i)),y*cos(Omega_k)*sin(i),
		0,sin(i),0,y*cos(i) };
	double X_dot_right[4] = { x_dot, y_dot, Omega_k_dot, I_dot };
	double X_dot[3] = { 0 };
	MatrixHelper m;
	if (!m.Matrix_Mult(R, 4, 3, X_dot_right, 1, 4, X_dot)) {
		cerr << "fault in GPS Matrix Mult" << endl;
		return 0;
	}
	sat.vel[0] = X_dot[0];
	sat.vel[1] = X_dot[1];
	sat.vel[2] = X_dot[2];
	//钟差钟漂计算(包含相对论)
	sat.clockbias = getClockbias(t) + getRelaCorr(E);
	sat.clockdrift = getClockdrift(t) + getRelaCorr_dot(E, n);
	//cout << this->satid << ":\n" << sat.x << endl;
	//cout << this->satid << ":\n" << sat.v << endl;
	//cout << Sigt << endl;
	//cout << this->satid << ":  " << u << '\t' << r << '\t' << i << endl;
	return 1;
}

int BDSEphemeris::getSatXvt(const GPSweek& t, PVT& sat)
{
	getBDST();
	//判断星历是否过期（1h）
	if (!isValid(t))
		return 0;
	double A = rootA * rootA;
	double tk = (t - toe).toSec();						//运行时间
	//toe时刻卫星平均角速度
	double n0 = sqrt(GM_BDS / (A*A*A));
	//观测时刻卫星平均角速度
	double n = n0 + deltaN;
	//平近点角
	double M = M0 + n * tk;
	//M = fmod(M + 2 * M_PI, 2 * M_PI);
	//偏近点角
	double E = M;
	if (!computeEk(M, E)) {
		return 0;
	}
	//真近点角
	double f = atan2(sqrt(1 - ecc * ecc)*sin(E), cos(E) - ecc);

	//纬度幅角
	double phi = f + omega;
	phi = fmod(phi, 2 * M_PI);
	//计算改正后的三个量
	double u = phi + cuc * cos(2 * phi) + cus * sin(2 * phi);
	double r = A * (1 - ecc * cos(E)) + crc * cos(2 * phi) + crs * sin(2 * phi);
	double i = i0 + i_dot * tk + cic * cos(2 * phi) + cis * sin(2 * phi);
	//计算卫星在轨道平面的坐标
	double pos = r * cos(u);
	double y = r * sin(u);

	//卫星速度计算
	double E_dot = n / (1 - ecc * cos(E));
	double phi_dot = sqrt((1 + ecc) / (1 - ecc))*
		pow((cos(f / 2) / cos(E / 2)), 2)*E_dot;

	double u_dot = 2 * (cus*cos(2 * phi) - cuc * sin(2 * phi))*phi_dot + phi_dot;
	double r_dot = A * ecc*sin(E)*E_dot + 2 * (crs*cos(2 * phi) - crc * sin(2 * phi))*phi_dot;
	double I_dot = i_dot + 2 * (cis*cos(2 * phi) - cic * sin(2 * phi))*phi_dot;
	double x_dot = r_dot * cos(u) - r * u_dot*sin(u);
	double y_dot = r_dot * sin(u) + r * u_dot*cos(u);

	MatrixHelper m;
	//若为MEO/IGSO轨道卫星
	if (satid.svprn >= 6 && satid.svprn <= 58) {
		//位置计算
		//升交点经度
		double Omega_k = Omega0 + (Omega_dot - Omega_e_BDS)*tk - Omega_e_BDS * bdstoe.weeksec;
		//卫星位置
		sat.pos.thearray[0] = pos * cos(Omega_k) - y * cos(i)*sin(Omega_k);
		sat.pos.thearray[1] = pos * sin(Omega_k) + y * cos(i)*cos(Omega_k);
		sat.pos.thearray[2] = y * sin(i);

		//速度计算
		double Omega_k_dot = Omega_dot - Omega_e_BDS;
		double X_dot[3] = { };
		//对位置求导
		double R[12] =
		{ cos(Omega_k),-sin(Omega_k)*cos(i),-(pos*sin(Omega_k) + y * cos(Omega_k)*cos(i)),y*sin(Omega_k)*sin(i),
			sin(Omega_k),cos(Omega_k)*cos(i),pos*cos(Omega_k) - y * sin(Omega_k)*cos(i),y*cos(Omega_k)*sin(i),
		0,sin(i),0,y*cos(i) };
		double X_dot_right[4] = { x_dot, y_dot, Omega_k_dot, I_dot };

		if (!m.Matrix_Mult(R, 4, 3, X_dot_right, 1, 4, X_dot)) {
			cerr << "fault in GPS Matrix Mult" << endl;
			return 0;
		}

		sat.vel.thearray[0] = X_dot[0];
		sat.vel.thearray[1] = X_dot[1];
		sat.vel.thearray[2] = X_dot[2];

	}
	//若为GEO轨道卫星
	else if (satid.svprn <= 5 || satid.svprn >= 59) {
		double Omega_k = Omega0 + Omega_dot * tk - Omega_e_BDS * bdstoe.weeksec;
		double X_GK[3] = { 0 };
		double X_K[3] = { 0 };
		double temp[9] = { 0 };
		double phi_t1 = Omega_e_BDS * tk;
		double phi_t2 = -5.0 / 180 * M_PI;

		X_GK[0] = pos * cos(Omega_k) - y * cos(i)*sin(Omega_k);
		X_GK[1] = pos * sin(Omega_k) + y * cos(i)*cos(Omega_k);
		X_GK[2] = y * sin(i);

		double Rz[9] = { cos(phi_t1),sin(phi_t1),0,
			-sin(phi_t1),cos(phi_t1),0,
			0,0,1 };

		double Rx[9] = { 1,0,0,
		0,cos(phi_t2),sin(phi_t2),
		0,-sin(phi_t2),cos(phi_t2) };
		if (!m.Matrix_Mult(Rz, 3, 3, Rx, 3, 3, temp)) {
			cerr << "fault in BDS Matrix Mult" << endl;
			return 0;
		}
		if (!m.Matrix_Mult(temp, 3, 3, X_GK, 1, 3, X_K)) {
			cerr << "fault in BDS Matrix Mult" << endl;
			return 0;
		}
		sat.pos.thearray[0] = X_K[0];
		sat.pos.thearray[1] = X_K[1];
		sat.pos.thearray[2] = X_K[2];

		//卫星速度计算
		double Omega_k_dot = Omega_dot;
		double X_GK_dot[3] = { 0 };
		double X_K_dot[3] = { 0 };
		double temp_dot[9] = { 0 };
		double result_dot1[3] = { 0 };
		double result_dot2[3] = { 0 };


		//对位置求导
		double R[12] =
		{ cos(Omega_k),-sin(Omega_k)*cos(i),-(pos*sin(Omega_k) + y * cos(Omega_k)*cos(i)),y*sin(Omega_k)*sin(i),
			sin(Omega_k),cos(Omega_k)*cos(i),pos*cos(Omega_k) - y * sin(Omega_k)*cos(i),y*cos(Omega_k)*sin(i),
		0,sin(i),0,y*cos(i) };
		double X_dot_right[4] = { x_dot, y_dot, Omega_k_dot, I_dot };

		if (!m.Matrix_Mult(R, 4, 3, X_dot_right, 1, 4, X_GK_dot)) {
			cerr << "fault in GPS Matrix Mult" << endl;
			return 0;
		}

		//对旋转矩阵求导
		double Rz_dot[9] = { -sin(phi_t1)*Omega_e_BDS,cos(phi_t1)*Omega_e_BDS,0,
		-cos(phi_t1)*Omega_e_BDS,-sin(phi_t1)*Omega_e_BDS,0,
		0,0,0 };
		for (short j = 0; j < 3; j++) {
			X_GK[j] = X_GK[j];
		}
		if (!m.Matrix_Mult(Rz_dot, 3, 3, Rx, 3, 3, temp_dot)) {
			cerr << "fault in BDS velocity Matrix" << endl;
		}
		if (!m.Matrix_Mult(temp_dot, 3, 3, X_GK, 1, 3, result_dot1)) {
			cerr << "fault in BDS velocity Matrix" << endl;
		}
		if (!m.Matrix_Mult(Rz, 3, 3, Rx, 3, 3, temp_dot)) {
			cerr << "fault in BDS velocity Matrix" << endl;
		}
		if (!m.Matrix_Mult(temp_dot, 3, 3, X_GK_dot, 1, 3, result_dot2)) {
			cerr << "fault in BDS velocity Matrix" << endl;
		}

		sat.vel.thearray[0] = result_dot1[0] + result_dot2[0];
		sat.vel.thearray[1] = result_dot1[1] + result_dot2[1];
		sat.vel.thearray[2] = result_dot1[2] + result_dot2[2];

	}
	else {
		cerr << "fault in BDS prn,PRN:" << this->satid.svprn << endl;
		return 0;
	}
	//钟差钟漂计算(包含相对论)
	sat.clockbias = getClockbias(t) + getRelaCorr(E);
	sat.clockdrift = getClockdrift(t) + getRelaCorr_dot(E, n);
	//cout << this->satid << ":\n" << sat.x << endl;
	//cout << this->satid << ":\n" << sat.v << endl;
	return 1;
}

GPSweek Ephemeris::getShootime(const GPSweek& trec, double psr)
{
	//信号发射时间计算
	short no_iter = 10;
	GPSweek traw, tsig(trec), temp;
	double dt = 0, clockbias = 0;
	for (short i = 0; i < no_iter; i++)
	{
		temp = tsig;
		traw = trec - psr / V_light;
		clockbias = getClockbias(traw);		//卫星钟差改正
		tsig = traw - clockbias;
		if (abs((temp - tsig).toSec()) < 1e-12)
		{
			//cout << clockbias << endl;
			return tsig;
		}
	}
	return tsig;
}

double GPSEphemeris::getRelaCorr(double Ek)
{
	double F = -4.442807633e-10;	//unit sec/sqrt(m)
	double tr = F * ecc * rootA * sin(Ek);
	return tr;
}

double BDSEphemeris::getRelaCorr(double Ek)
{
	double F = -4.442807633e-10;	//unit sec/sqrt(m)
	double tr = F * ecc * rootA * sin(Ek);
	return tr;
}

double GPSEphemeris::getRelaCorr_dot(double Ek, double n) {
	double F = -4.442807633e-10;	//unit sec/sqrt(m)
	double E_dot = n / (1 - ecc * cos(Ek));
	double tr = F * ecc * rootA * sin(Ek) * E_dot;
	return tr;
}

double BDSEphemeris::getRelaCorr_dot(double Ek, double n) {
	double F = -4.442807633e-10;	//unit sec/sqrt(m)
	double E_dot = n / (1 - ecc * cos(Ek));
	double tr = F * ecc * rootA * sin(Ek) * E_dot;
	return tr;
}

bool Ephemeris::computeEk(double M, double &E)
{
	//偏近点角
	E = M;
	short no_iter = 10;
	for (int j = 0; j < no_iter; ++j) {
		double E_old = E;
		E = M + ecc * sin(E);
		double dE = fmod(E - E_old, 2 * M_PI);
		if (abs(dE) <= 1e-12)break;

		if (j == no_iter) {
			cerr << "fault in BDS sat position" << endl;
			return false;
		}
	}
	E = fmod(E + 2 * M_PI, 2 * M_PI);
	return true;
}

double Ephemeris::getClockbias(const GPSweek& t)
{
	//计算卫星钟差
	double dt = (t - toc).toSec();
	double tcorr = (af2 * dt + af1)*dt + af0;
	return tcorr;
}

double Ephemeris::getClockdrift(const GPSweek& t)
{
	return(af1 + 2 * af2 * (t - toc).toSec());
}

Triple GPSEphemeris::getEarthRotCorr(double traveltime, const Triple& Xraw)
{
	MatrixHelper m;
	double X[3] = { Xraw[0], Xraw[1], Xraw[2] };
	double X_rot[3] = { 0 };
	double rotation = Omega_e_GPS * traveltime;

	double Rotate[9] = { cos(rotation), sin(rotation), 0,
	-sin(rotation), cos(rotation), 0,
	0, 0, 1 };

	m.Matrix_Mult(Rotate, 3, 3, X, 1, 3, X_rot);
	Triple result(X_rot[0], X_rot[1], X_rot[2]);
	return result;
}

Triple BDSEphemeris::getEarthRotCorr(double traveltime, const Triple& Xraw)
{
	MatrixHelper m;
	double X[3] = { Xraw[0], Xraw[1], Xraw[2] };
	double X_rot[3] = { 0 };
	double rotation = Omega_e_BDS * traveltime;

	double Rotate[9] = { cos(rotation), sin(rotation), 0,
	-sin(rotation), cos(rotation), 0,
	0, 0, 1 };

	m.Matrix_Mult(Rotate, 3, 3, X, 1, 3, X_rot);
	Triple result(X_rot[0], X_rot[1], X_rot[2]);
	return result;
}

bool GPSEphemeris::isValid(const GPSweek& t)const {
	if (abs((t - this->toe).toSec()) <= 7200)
		return true;
	else
		return false;
}
bool BDSEphemeris::isValid(const GPSweek& t)const {
	if (abs((t - this->toe).toSec()) <= 3600)
		return true;
	else
		return false;
}
void BDSEphemeris::getBDST() {
	TimeTrans::GPST2BDST(toe, bdstoe);
}
shared_ptr<Ephemeris> GPSEphemeris::clone()const {
	shared_ptr<Ephemeris> p(new GPSEphemeris(*this));
	return p;
}
shared_ptr<Ephemeris> BDSEphemeris::clone()const {
	shared_ptr<Ephemeris> p(new BDSEphemeris(*this));
	return p;
}