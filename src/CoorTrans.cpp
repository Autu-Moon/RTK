#include"CoorTrans.h"
//笛卡尔坐标到大地坐标
//笛卡尔坐标转大地坐标
//输出B L 为弧度格式
bool CoorTrans::XYZToBLH(const Triple &xyz, Triple &blh)
{
	Triple blh_t;
	double delta_Z = ref.esq * xyz[2];
	blh_t[1] = atan2(xyz[1], xyz[0]);
	short no_iteration = 15;//设置最大迭代次数
	double temp = 0;
	for (int i = 0; i < no_iteration; i++)
	{
		temp = delta_Z;
		blh_t[0] = atan2((xyz[2] + delta_Z), sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]));
		double sinB = (xyz[2] + delta_Z) / sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + pow(xyz[2] + delta_Z, 2));
		double N = ref.a / sqrt(1 - ref.esq*sinB*sinB);
		blh_t[2] = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + pow((xyz[2] + delta_Z), 2)) - N;
		delta_Z = N * ref.esq * sinB;
		if (abs(temp - delta_Z) < 1e-10) {
			blh = blh_t;
			return true;
			break;
		}
	}
#ifdef DEBUG
	std::cerr << "XYZ to BLH fault" << std::endl;
#endif
	return false;
}
//大地坐标到笛卡尔坐标
//输入B L 为弧度格式
void CoorTrans::BLHToXYZ(const Triple& blh, Triple &xyz)
{
	double N = ref.computeR_N(blh[0]);
	xyz[0] = (N + blh[2]) * cos(blh[0]) * cos(blh[1]);
	xyz[1] = (N + blh[2]) * cos(blh[0]) * sin(blh[1]);
	xyz[2] = (N * (1 - ref.esq) + blh[2]) * sin(blh[0]);
}
bool CoorTrans::XYZToNEU(const Triple &xyzSat, const Triple &xyzRec, Triple &NEU, double &elev)
{
	MatrixHelper m;
	Triple blhRec;
	/*if (XYZToBLH(xyzSat, blhRec) == false)
		return false;*/
	XYZToBLH(xyzRec, blhRec);
	double Neu[3] = { 0 };
	double Trans[9] = { -sin(blhRec[0])*cos(blhRec[1]),-sin(blhRec[0])*sin(blhRec[1]),cos(blhRec[0]),
	-sin(blhRec[1]),cos(blhRec[1]),0,
	cos(blhRec[0])*cos(blhRec[1]),cos(blhRec[0])*sin(blhRec[1]),sin(blhRec[0]) };
	double Xx[3] = { xyzSat[0] - xyzRec[0],xyzSat[1] - xyzRec[1], xyzSat[2] - xyzRec[2] };

	if (!m.Matrix_Mult(Trans, 3, 3, Xx, 1, 3, Neu))
		return false;
	double hos_dis = sqrt(Neu[0] * Neu[0] + Neu[1] * Neu[1]);
	elev = atan(Neu[2] / hos_dis);
	NEU[0] = Neu[0];
	NEU[1] = Neu[1];
	NEU[2] = Neu[2];

	return true;
}
bool CoorTrans::XYZToNEU(const Triple &xyzSat, const Triple &xyzRec, Triple &NEU)
{
	MatrixHelper m;
	Triple blhRec;
	/*if (XYZToBLH(xyzSat, blhRec) == false)
		return false;*/
	XYZToBLH(xyzRec, blhRec);
	double Neu[3] = { 0 };
	double Trans[9] = { -sin(blhRec[0])*cos(blhRec[1]),-sin(blhRec[0])*sin(blhRec[1]),cos(blhRec[0]),
	-sin(blhRec[1]),cos(blhRec[1]),0,
	cos(blhRec[0])*cos(blhRec[1]),cos(blhRec[0])*sin(blhRec[1]),sin(blhRec[0]) };
	double Xx[3] = { xyzSat[0] - xyzRec[0],xyzSat[1] - xyzRec[1], xyzSat[2] - xyzRec[2] };

	if (!m.Matrix_Mult(Trans, 3, 3, Xx, 1, 3, Neu))
		return false;

	NEU[0] = Neu[0];
	NEU[1] = Neu[1];
	NEU[2] = Neu[2];

	return true;
}
bool CoorTrans::XYZToNEU(const Triple &xyzSat, const Triple &xyzRec, double &elev)
{
	MatrixHelper m;
	Triple blhRec;
	/*if (XYZToBLH(xyzSat, blhRec) == false)
		return false;*/
	XYZToBLH(xyzRec, blhRec);
	double Neu[3] = { 0 };
	double Trans[9] = { -sin(blhRec[0])*cos(blhRec[1]),-sin(blhRec[0])*sin(blhRec[1]),cos(blhRec[0]),
	-sin(blhRec[1]),cos(blhRec[1]),0,
	cos(blhRec[0])*cos(blhRec[1]),cos(blhRec[0])*sin(blhRec[1]),sin(blhRec[0]) };
	double Xx[3] = { xyzSat[0] - xyzRec[0],xyzSat[1] - xyzRec[1], xyzSat[2] - xyzRec[2] };

	if (!m.Matrix_Mult(Trans, 3, 3, Xx, 1, 3, Neu)) {
		elev = 0;
		return false;
	}

	double hos_dis = sqrt(Neu[0] * Neu[0] + Neu[1] * Neu[1]);
	elev = atan(Neu[2] / hos_dis);

	return true;
}
