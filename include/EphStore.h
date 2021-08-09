#pragma once
#include"Ephemeris.h"
#include<map>
using std::map;
//储存所有卫星星历 并提供getXvt被调用
class EphStore
{
	//储存GPS+BDS星历
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
	* @Description: 实现卫星位置速度，信号发射时间，地球自转改正的计算
	* @Param: 卫星号 伪距观测值 接收机表面时
	* @return: 卫星 位置 速度 钟差 钟漂
	*/
	int getXvt(const SatID &sat, double psr, const GPSweek& t, PVT &xvt)const;

	/*
	* @Description: 实现卫星位置速度，信号发射时间，地球自转改正的计算
	* @Param: 卫星号 观测值数据
	* @return: 卫星 位置 速度 钟差 钟漂
	*/
	int getXvt(const SatID &sat, ObsStore &obs, PVT &xvt)const;

	/*
	* @Description: 实现卫星位置速度，信号发射时间，地球自转改正的计算
	* @Param: 卫星号 信号发射时间
	* @return: 卫星 位置 速度 钟差 钟漂
	*/
	int getXvt(const SatID &sat, GPSweek t, PVT &xvt)const;
};