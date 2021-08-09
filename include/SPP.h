#pragma once
#include"EphStore.h"
#include"ObsStore.h"
#include"UserObsStore.h"
#include"TropoCorr.h"
#include<fstream>
#include<iomanip>
using namespace std;
//#define DEBUG_SPP 1
//单点定位类 储存接收机/卫星PVT 并进行单点定位计算
class SPP 
{
public:
	//-------------field-------------
	Triple		realRecPos = Triple(-2267804.5263,5009342.3723,3220991.8632);		//位置真值
	ObsID		SPPobsID;															//观测值类型
	/*
	 * @Description: no data struct, UserObsStore provide data needed by SPP directly
	 */
	struct SPPresult
	{
		GPSweek		epoch;
		Triple		PositionXYZ;														//unit:		m	m	m 
		Triple		deltaNEU;														//unit:		m	m	m 
		Triple		PositionBLH;														//unit:		deg	deg	deg
		double		PositionSigma;
		double		PDOP;

		Triple		VelocityXYZ;
		double		VelocitySigma;

		int			SYSnum[SYScount];
		void Empty() {
			epoch = GPSweek(0, 0);
			PositionXYZ = PositionBLH = deltaNEU = VelocityXYZ = Triple(0, 0, 0);
			PositionSigma = 0;
			PDOP = 0;
			VelocitySigma = 0;
			std::fill(SYSnum, SYSnum + SYScount, 0);
		}
	};
	SPPresult	result;															//储存单点定位结果
	//---------set function---------
	void		setObsID(const ObsID &id);										//接受观测值类型
	//--------solve function--------
	/*
	* @Description: SPP
	* @Param: UserObsStore
	*/
	//单点定位（最小二乘）
	int			solvePosition_LS(const UserObsStore& usrStore);
	int			modelPosSingleSys_LS(const UserObsStore& usrStore, SYSTEM sys);		//单系统
	int			modelPosDoubleSys_LS(const UserObsStore& usrStore);					//双系统
	//单点测速（最小二乘）
	int			solveVel_LS(const UserObsStore& usrStore);		
	int			modelVelocity_LS(const UserObsStore& usrStore);						//双系统
	//总调用函数
	int			Solve(const UserObsStore& usrStore);
	int			Solve(const UserObsStore& usrStore,ostream& os);
	//-----------function-----------
	/*
	* @Description: << result to ostream
	* @Param: ostream
	* @return: void
	*/
	void		dumpResult(ostream& os)const;
	/*
	 * @return: position in xyz
	 */
	Triple		getPos()const;

	void		setRealPos(const Triple& pos);
};
