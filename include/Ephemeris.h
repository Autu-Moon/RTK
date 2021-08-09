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


//储存单个卫星的导航电文值
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
	double ecc;		//离心率

	double omega;	//近地点幅角
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
	* @Description: 计算卫星信号发射时间
	* @Param: 接收机收到卫星信号时间，伪距测量值
	* @return: 卫星信号发射时间
	*/
	virtual GPSweek getShootime(const GPSweek& trec, double psr);

	/*
	* @Description: 计算卫星钟差
	* @Param: 卫星信号发射时间
	* @return: 卫星钟差
	*/
	virtual double getClockbias(const GPSweek& t);

	/*
	* @Description: 计算卫星钟漂
	* @Param: 信号发射时间
	* @return: 卫星钟漂
	*/
	virtual double getClockdrift(const GPSweek& t);

	/*
	* @Description: 计算相对论改正
	* @Param: 偏近点角
	* @return: 相对论改正值
	*/
	virtual double getRelaCorr(double Ek) = 0;

	/*
	* @Description: 计算相对论改正的导数
	* @Param: 偏近点角 观测时刻卫星平均角速度
	* @return: 相对论改正的导数
	*/
	virtual double getRelaCorr_dot(double Ek, double n) = 0;

	/*
	* @Description: 计算卫星的位置 速度 钟差 钟漂（包含相对论改正值）
	* @Param: 信号发射时间 卫星PVT结构体
	* @return: 0 for fault, 1 for ok
	*/
	virtual int getSatXvt(const GPSweek& t, PVT& sat) = 0;

	/*
	* @Description: 对卫星位置进行地球自转改正
	* @Param: 信号传播时间 未经改正的卫星位置
	* @return: 改正后的卫星位置
	*/
	virtual Triple getEarthRotCorr(double traveltime, const Triple &Xraw) = 0;

	//-----------function------------
	/*
	* @Description: 计算偏近点角
	* @Param: 平近点角 偏近点角
	* @return: 是否成功计算
	*/
	virtual bool computeEk(double M, double &E);
	/*
	* @Description: 返回智能指针
	* @Param:
	* @return: 指向该类的智能指针
	*/
	virtual shared_ptr<Ephemeris> clone() const = 0;
	/*
	* @Description: judge if this eph is valid for t
	* @Param: gpst
	* @return: bool
	*/
	virtual bool isValid(const GPSweek& t) const = 0;
};

//GPS星历类 并实现卫星位置速度计算和地球自转改正
class GPSEphemeris :public Ephemeris
{
public:
	//------------field-------------
	/*GPS Tgd[0]*/
	unsigned long IODE;
	//---------get function---------
	int getSatXvt(const GPSweek& t, PVT& sat)override;
	Triple getEarthRotCorr(double traveltime, const Triple &Xraw)override;
	double getRelaCorr(double Ek)override;					//计算相对论改正
	double getRelaCorr_dot(double Ek, double n)override;	//计算相对论改正的导数
	//-----------function-----------
	bool isValid(const GPSweek& t)const override;
	shared_ptr<Ephemeris> clone()const override;
	//----------costructor----------
	GPSEphemeris() = default;
	GPSEphemeris(SYSTEM s, short prn) :Ephemeris(s, prn) {};
};

//BDS星历结构体 并实现卫星位置速度计算和地球自转改正
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
	double getRelaCorr(double Ek)override;					//计算相对论改正
	double getRelaCorr_dot(double Ek, double n)override;	//计算相对论改正的导数
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


