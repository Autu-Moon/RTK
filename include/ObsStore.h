#pragma once
#include"TimeTrans.h"
#include"SatID.h"
#include"PVT.h"
#include"Ephemeris.h"
#include<map>
#include<algorithm>
#include<deque>
#include<iomanip>
#include"SatNavBasic.h"

class ObsID
{
public:
	std::string Marker;					//测站名
	double interval;
	std::deque<FreqBand> GPSfremark;
	std::deque<FreqBand> BDSfremark;
	ObsID() {};
	void loadGPS(FreqBand g) {
		GPSfremark.push_back(g);
		//接收后排序
		sort(GPSfremark.begin(), GPSfremark.end());
	}
	void loadBDS(FreqBand b) {
		BDSfremark.push_back(b);
		//接收后排序
		sort(BDSfremark.begin(), BDSfremark.end());
	}
	ObsID& operator=(const ObsID& right) {
		this->GPSfremark.clear();
		this->BDSfremark.clear();
		for (auto g : right.GPSfremark) {
			this->GPSfremark.push_back(g);
		}
		for (auto b : right.BDSfremark) {
			this->BDSfremark.push_back(b);
		}
		return *this;
	}
};

//储存单历元下所有卫星的观测值
struct ObsData
{
	SatID satid;								//卫星号
	//L1   L2    L5
	//B1I  B3I
	double Psr[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };		//伪距
	double Adr[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };		//载波相位
	double Dop[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };		//多普勒
	double CNR[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };		//载波信噪比
	//double Psr_sigma[3];						//伪距误差
	//double Adr_sigma[3];						//载波相位误差
	ObsData() {};
	ObsData(SYSTEM s, short prn) :satid(s, prn) {};
	void reset() {
		std::fill(Psr, Psr + OBSERVATIONLENGTH, OBSEMPTY);
		std::fill(Adr, Adr + OBSERVATIONLENGTH, OBSEMPTY);
		std::fill(Dop, Dop + OBSERVATIONLENGTH, OBSEMPTY);
		std::fill(CNR, CNR + OBSERVATIONLENGTH, OBSEMPTY);
		satid.reset();
	}
	static bool isEmpty(const double& obs) {
		if (abs(obs - OBSEMPTY) < 1e-6) {
			return true;
		}
		else
			return false;
	}
};
class ObsStore
{
public:
	//------------field-------------
	//储存单历元下单个卫星的双频观测值

	typedef std::map<SatID, ObsData> EpochObsMap;
	////储存所有卫星观测值的map
	//typedef std::map<GPSweek, std::map<SatID, ObsData>> EpochObsMap;
	//
	//EpochObsMap EpobsMap;
	EpochObsMap EpobsMap;
	ObsID obsid;
	GPSweek epoch;
	//-----------function-----------
	/*
	* @Description:  add data from a sat to the map
	* @Param: SatID,ObsData
	* @return: void
	*/
	void addObs(const SatID& sat,const ObsData& obsdata);
	/*
	* @Description: cout all of the map to Console
	*/
	void dump(std::ostream& os)const;
	/*
	* @Description: clear the map
	*/
	void clear();
	//----------costructor----------
	ObsStore() = default;
};
