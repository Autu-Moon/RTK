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
	std::string Marker;					//��վ��
	double interval;
	std::deque<FreqBand> GPSfremark;
	std::deque<FreqBand> BDSfremark;
	ObsID() {};
	void loadGPS(FreqBand g) {
		GPSfremark.push_back(g);
		//���պ�����
		sort(GPSfremark.begin(), GPSfremark.end());
	}
	void loadBDS(FreqBand b) {
		BDSfremark.push_back(b);
		//���պ�����
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

//���浥��Ԫ���������ǵĹ۲�ֵ
struct ObsData
{
	SatID satid;								//���Ǻ�
	//L1   L2    L5
	//B1I  B3I
	double Psr[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };		//α��
	double Adr[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };		//�ز���λ
	double Dop[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };		//������
	double CNR[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };		//�ز������
	//double Psr_sigma[3];						//α�����
	//double Adr_sigma[3];						//�ز���λ���
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
	//���浥��Ԫ�µ������ǵ�˫Ƶ�۲�ֵ

	typedef std::map<SatID, ObsData> EpochObsMap;
	////�����������ǹ۲�ֵ��map
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
