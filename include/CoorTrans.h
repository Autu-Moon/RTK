#pragma once
#include"MatrixHelper.h"
#include"PVT.h"
#include"RefFrame.h"
class CoorTrans
{
	//用于同一坐标系统下的坐标形式转换
private:
	RefFrame ref;
public:
	CoorTrans(RefFrame::RefFrameType reftype=RefFrame::WGS84):ref(reftype){}

	bool XYZToBLH(const Triple &xyz, Triple &blh);
	void BLHToXYZ(const Triple &blh, Triple &xyz);
	bool XYZToNEU(const Triple &xyzSat, const Triple &xyzRec, Triple &NEU, double &elev);
	bool XYZToNEU(const Triple &xyzSat, const Triple &xyzRec, Triple &NEU);
	bool XYZToNEU(const Triple &xyzSat, const Triple &xyzRec, double &elev);
};
