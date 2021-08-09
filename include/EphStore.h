#pragma once
#include"Ephemeris.h"
#include<map>
using std::map;
//���������������� ���ṩgetXvt������
class EphStore
{
	//����GPS+BDS����
public:
	using TimeEphMap = std::map<GPSweek, shared_ptr<Ephemeris>>;
	using EphMap = std::map<SatID, TimeEphMap>;
	EphMap ephMap;
	/*
	* @Description: add eph to the main map
	* @Param: eph to add
	*/
	virtual void addEphemeris(const Ephemeris& eph);

	/*
	* @Description: find the nearest ephemeris
	* @Param: satid gpst
	* @return: shared_ptr to the ephemeris
	*/
	const shared_ptr<Ephemeris> findEph(const SatID& sat, const GPSweek& t)const;

	/*
	* @Description: get the TimeEphMap of sat
	* @Param: SatID
	* @return: bool ; TimeEphMap
	*/
	bool getTimeEphMap(const SatID& sat, TimeEphMap& tem)const;

	/*
	* @Description: ʵ������λ���ٶȣ��źŷ���ʱ�䣬������ת�����ļ���
	* @Param: ���Ǻ� α��۲�ֵ ���ջ�����ʱ
	* @return: ���� λ�� �ٶ� �Ӳ� ��Ư
	*/
	int getXvt(const SatID &sat, double psr, const GPSweek& t, PVT &xvt)const;

	/*
	* @Description: ʵ������λ���ٶȣ��źŷ���ʱ�䣬������ת�����ļ���
	* @Param: ���Ǻ� �۲�ֵ����
	* @return: ���� λ�� �ٶ� �Ӳ� ��Ư
	*/
	int getXvt(const SatID &sat, ObsStore &obs, PVT &xvt)const;

	/*
	* @Description: ʵ������λ���ٶȣ��źŷ���ʱ�䣬������ת�����ļ���
	* @Param: ���Ǻ� �źŷ���ʱ��
	* @return: ���� λ�� �ٶ� �Ӳ� ��Ư
	*/
	int getXvt(const SatID &sat, GPSweek t, PVT &xvt)const;
};