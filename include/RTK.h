#pragma once
#include"MathUtils.h"
#include"SatNavBasic.h"
#include"EphStore.h"
#include"ObsStore.h"
#include"UserObsStore.h"
#include"MyMatrix.h"
#include"LeastSquare.h"
#include"lambda.h"
#include<iostream>
#include<fstream>
#include<iomanip>
//#define DEBUG_RTK 1
class RTK
{
public:
	/*
	 * @Description: realize RTK 
	 * @Param: pivot: base satellite
	 * @Param: norm: none base satellite
	 * @Param: base: base station
	 * @Param: rover: rover station
	 * @Param: doudiff: (nr - nb) - (pr - pb)
	 */
	//------------field-------------
	Triple				realRecPos;										//基准站坐标
	ObsID				RTKobsID;										//观测值类型
	static double		RatioThreshold;

	SatID				PivotSat[SYScount];								//基准星 GPS BDS
	static const Triple	BaseRec;										//基准站已知坐标
	/*
	 * @Description: provide data needed by RTK
	 * obscode(_C = Psr, _L = Adr)
	 * SYSTEM
	 * FreqBand
	 * SatID - double 
	 */
	using				DouDiffMap = map<ObsCode, map<SYSTEM, map<FreqBand, map<SatID, double>>>>;
	DouDiffMap			doudiffMap;
	struct RTKresult
	{
		GPSweek		epoch;
		Triple		PositionXYZ;									//unit:		m	m	m 
		Triple		deltaNEU;										//unit:		m	m	m 
		Triple		PositionBLH;									//unit:		deg	deg	deg
		Triple		Baseline;
		bool		fixed;
		double		PositionSigma;
		double		PDOP;

		int			SatNum[SYScount];
		void Empty() {
			epoch = GPSweek(0, 0);
			PositionXYZ = PositionBLH = deltaNEU = Baseline = Triple(0, 0, 0);
			PositionSigma = 0;
			PDOP = 0;
			SatNum[GPS] = 0;
			SatNum[BDS] = 0;
			fixed = false;
		}
	};
	RTKresult	result;													//储存单点定位结果
	//---------set function---------
	void		setObsID(const ObsID &id);								//接受观测值类型
	/*
	 * @Description: set Double difference observation
	 * @Param: base UserObsStore, rover UserObsStore, system
	 * @return: Triple[0] for number of psr, Triple[1] for number of adr
	 */
	Triple		setDouDiffMap(UserObsStore& base, UserObsStore& rover, SYSTEM sys);
	/*
	 * @Description: set realRecPos
	 * @Param: Triple
	 * @return: void
	 */
	void		setRealRecPos(const Triple& real);
	//--------solve function--------
	/*
	* @Description: RTK
	* @Param: UserObsStore
	*/
	//RTK（最小二乘）
	/*
	 * @Description: model and solve RTK
	 * @Param: base UserObsStore, rover UserObsStore, 
	 * @Param: ObsCount[0] for  number of psr, [1] for number of adr
	 * @Param: coordinate of rover
	 * @return: 1 for success
	 */
	int				modelPosition_LS(const UserObsStore& base, const UserObsStore& rover, const Triple& ObsCount, bool fixed,
						const Triple& RoverX, const Triple& BaseX = BaseRec);
	/*
	 * @Description: model weight Matrix by CNR or elev
	 * @Param: 
	 * @return: weight Matrix 
	 */
	mym::Matrixd	modelCovDx(const UserObsStore& base, const UserObsStore& rover, 
						const ObsCode& code,const SYSTEM& sys,const FreqBand& freq);
	//总调用函数
	/*
	 * @Description: main function to process RTK
	 * @Param: base UserObsStore, rover UserObsStore, coordinate of rover,(stream)
	 */
	int				Solve(UserObsStore& base, UserObsStore& rover,bool fixed, const Triple& RoverX,  const Triple& BaseX = BaseRec);
	int				Solve(UserObsStore& base, UserObsStore& rover, ostream& os,bool fixed, const Triple& RoverX,  const Triple& BaseX = BaseRec);
	//-----------function-----------
	/*
	 * @Description: 
	 * @Param: 
	 * @return: 
	 */
	void		match();
	/*
	 * @Description: select Pivot satellite
	 * @Param: base rover
	 * @return:	base SatID 
	 */
	bool		selectPivotSat(UserObsStore& base, UserObsStore& rover, SYSTEM sys, SatID& satid);
	/*
	 * @Description: compare obs of two station, and erase rover sats that base doesn't have
	 * @Param: base rover
	 * @return: void
	 */
	void		selectRoverObs(UserObsStore& base, UserObsStore& rover);

	/*
	 * @Description: compare obs of two station, and erase base sats that rover doesn't have 
	 * @Param: base rover
	 * @return: void
	 */
	void		selectBaseObs(UserObsStore& base, UserObsStore& rover);
	/*
	 * @Description: compute double differece observation
	 *		(nr - nb) - (pr - pb)
	 * @Param: pivot-base pivot-rover norm-base norm-rover
	 * @return: double differece observation
	 */
	double		computeDouDiff(double pb, double pr, double nb, double nr);
	/*
	 * @Description: compute prioi sigma of a observation based on elevating angle
	 * @Param: elevating angle, ObsCode(_C _L)
	 * @return: prioi sigma
	 */
	double		weightElev(double elev,const ObsCode& code);
	/*
	 * @Description: fix ambiguity by lambda method
	 * @Param: covDx-Q, correct number(N)-dx, num of fixed solution-fixedNum, 
	 * @return: bool
	 */
	bool		fixAmbiguity(const mym::Matrixd& Q, const mym::Matrixd dx, int fixedNum, mym::Matrixd& N, double RatThre = RatioThreshold);
	/*
	* @Description: << result to ostream
	* @Param: ostream
	* @return: void
	*/
	void		dumpResult(std::ostream& os);
	/*
	* @Description: << double diff map to ostream
	* @Param: ostream
	* @return: void
	*/
	void		dumpDouDiff(std::ostream& os);
	Triple		computeCoeffB(const Triple& rec, const Triple& sat);
	//----------costructor----------
	RTK() = default;
	~RTK() = default;
};

