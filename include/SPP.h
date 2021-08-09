#pragma once
#include"EphStore.h"
#include"ObsStore.h"
#include"UserObsStore.h"
#include"TropoCorr.h"
#include<fstream>
#include<iomanip>
using namespace std;
//#define DEBUG_SPP 1
//���㶨λ�� ������ջ�/����PVT �����е��㶨λ����
class SPP 
{
public:
	//-------------field-------------
	Triple		realRecPos = Triple(-2267804.5263,5009342.3723,3220991.8632);		//λ����ֵ
	ObsID		SPPobsID;															//�۲�ֵ����
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
	SPPresult	result;															//���浥�㶨λ���
	//---------set function---------
	void		setObsID(const ObsID &id);										//���ܹ۲�ֵ����
	//--------solve function--------
	/*
	* @Description: SPP
	* @Param: UserObsStore
	*/
	//���㶨λ����С���ˣ�
	int			solvePosition_LS(const UserObsStore& usrStore);
	int			modelPosSingleSys_LS(const UserObsStore& usrStore, SYSTEM sys);		//��ϵͳ
	int			modelPosDoubleSys_LS(const UserObsStore& usrStore);					//˫ϵͳ
	//������٣���С���ˣ�
	int			solveVel_LS(const UserObsStore& usrStore);		
	int			modelVelocity_LS(const UserObsStore& usrStore);						//˫ϵͳ
	//�ܵ��ú���
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
