#pragma once
#include<iostream>
#include<vector>
#include"SatNavBasic.h"
#include"TimeTrans.h"
#include"CoorTrans.h"
#include"MathUtils.h"
#include"ObsStore.h"
#include"MatrixHelper.h"
#include"SatID.h"
using namespace std;
using std::vector;
const double GM_GPS = 3.986005e14;					//earth's universal gravitational parameter 
const double GM_BDS = 3.986004418e14;				//unit:m^3/s^2

const double Omega_e_GPS = 7.2921151467e-5;			// earth rotation rate, 
const double Omega_e_BDS = 7.2921150e-5;			//unit:rad/s


//���浥�����ǵĵ�������ֵ
class Ephemeris
{
public:
	//------------field-------------
	GPSweek toe;
	SatID satid;
	double tow;
	double rootA;		//Semi-major axis (m) 
	double deltaN;
	double M0;
	double ecc;		//������

	double omega;	//���ص����
	double cuc;
	double cus;
	double crc;
	double crs;
	double cic;
	double cis;
	double i0;															//unit:radian
	double i_dot;														//unit:radian/sec
	double Omega0;														//unit:radian
	double Omega_dot;													//unit:radian/sec

	GPSweek toc;
	double af0;
	double af1;
	double af2;
	short health;		//Sat health
	double Tgd[EPHTGDLENGTH];
	//----------costructor----------
	Ephemeris() = default;
	Ephemeris(SYSTEM s, short prn) :satid(s, prn) {};

	virtual ~Ephemeris() {};
	//---------get function----------

	/*
	* @Description: ���������źŷ���ʱ��
	* @Param: ���ջ��յ������ź�ʱ�䣬α�����ֵ
	* @return: �����źŷ���ʱ��
	*/
	virtual GPSweek getShootime(const GPSweek& trec, double psr);

	/*
	* @Description: ���������Ӳ�
	* @Param: �����źŷ���ʱ��
	* @return: �����Ӳ�
	*/
	virtual double getClockbias(const GPSweek& t);

	/*
	* @Description: ����������Ư
	* @Param: �źŷ���ʱ��
	* @return: ������Ư
	*/
	virtual double getClockdrift(const GPSweek& t);

	/*
	* @Description: ��������۸���
	* @Param: ƫ�����
	* @return: ����۸���ֵ
	*/
	virtual double getRelaCorr(double Ek) = 0;

	/*
	* @Description: ��������۸����ĵ���
	* @Param: ƫ����� �۲�ʱ������ƽ�����ٶ�
	* @return: ����۸����ĵ���
	*/
	virtual double getRelaCorr_dot(double Ek, double n) = 0;

	/*
	* @Description: �������ǵ�λ�� �ٶ� �Ӳ� ��Ư����������۸���ֵ��
	* @Param: �źŷ���ʱ�� ����PVT�ṹ��
	* @return: 0 for fault, 1 for ok
	*/
	virtual int getSatXvt(const GPSweek& t, PVT& sat) = 0;

	/*
	* @Description: ������λ�ý��е�����ת����
	* @Param: �źŴ���ʱ�� δ������������λ��
	* @return: �����������λ��
	*/
	virtual Triple getEarthRotCorr(double traveltime, const Triple &Xraw) = 0;

	//-----------function------------
	/*
	* @Description: ����ƫ�����
	* @Param: ƽ����� ƫ�����
	* @return: �Ƿ�ɹ�����
	*/
	virtual bool computeEk(double M, double &E);
	/*
	* @Description: ��������ָ��
	* @Param:
	* @return: ָ����������ָ��
	*/
	virtual shared_ptr<Ephemeris> clone() const = 0;
	/*
	* @Description: judge if this eph is valid for t
	* @Param: gpst
	* @return: bool
	*/
	virtual bool isValid(const GPSweek& t) const = 0;
};

//GPS������ ��ʵ������λ���ٶȼ���͵�����ת����
class GPSEphemeris :public Ephemeris
{
public:
	//------------field-------------
	/*GPS Tgd[0]*/
	unsigned long IODE;
	//---------get function---------
	int getSatXvt(const GPSweek& t, PVT& sat)override;
	Triple getEarthRotCorr(double traveltime, const Triple &Xraw)override;
	double getRelaCorr(double Ek)override;					//��������۸���
	double getRelaCorr_dot(double Ek, double n)override;	//��������۸����ĵ���
	//-----------function-----------
	bool isValid(const GPSweek& t)const override;
	shared_ptr<Ephemeris> clone()const override;
	//----------costructor----------
	GPSEphemeris() = default;
	GPSEphemeris(SYSTEM s, short prn) :Ephemeris(s, prn) {};
};

//BDS�����ṹ�� ��ʵ������λ���ٶȼ���͵�����ת����
class BDSEphemeris :public Ephemeris
{
public:
	//------------field-------------
	/*Tgd[0]=Tgd13 Tgd[1]=Tgd23;*/
	unsigned long AODE;

	GPSweek bdstoe;  //field"toe" is in BDST, "bdstoe" is in BDST
	//---------get function---------
	int getSatXvt(const GPSweek& t, PVT& sat)override;
	Triple getEarthRotCorr(double traveltime, const Triple &Xraw)override;
	double getRelaCorr(double Ek)override;					//��������۸���
	double getRelaCorr_dot(double Ek, double n)override;	//��������۸����ĵ���
	//-----------function-----------
	/*
	* @Description: trans bdst to gpst ,stored in field "gpstoe"
	*/
	void getBDST();
	bool isValid(const GPSweek& t)const override;
	shared_ptr<Ephemeris> clone()const override;
	//----------costructor----------
	BDSEphemeris() = default;
	BDSEphemeris(SYSTEM s, short prn) :Ephemeris(s, prn) {};
};


