#include"EphStore.h"
int EphStore::getXvt(const SatID &sat, ObsStore &obsStore, PVT &xvt) const{
	//根据给出的卫星号，从储存的星历中算出卫星位置速度
	//信号发射时间
	//0 for no such sat
	int result = 0;
	EphMap::const_iterator iteph = ephMap.find(sat);											//判断是否有该卫星星历
	if (iteph == ephMap.end())
		return result;
	ObsStore::EpochObsMap::iterator itobs = obsStore.EpobsMap.find(sat);
	if (itobs == obsStore.EpobsMap.end()) {
		return result;
	}
	result = getXvt(sat, itobs->second.Psr[0], obsStore.epoch, xvt);
	return result;
}
int EphStore::getXvt(const SatID &sat, double psr, const GPSweek& t, PVT &xvt)const{
	//根据给出的卫星号，从储存的星历中算出卫星位置速度
	//信号发射时间
	//0 for no such sat
	int result = 0;
	EphMap::const_iterator iteph = ephMap.find(sat);											//判断是否有该卫星星历
	if (iteph == ephMap.end())
		return result;
	//find proper Ephemeris
	shared_ptr<Ephemeris> eph = findEph(sat,t);
	if (!eph) {
#ifdef DEBUG
		std::cerr << sat << " no eph in " << t << std::endl;
#endif
		return result;
	}
	//信号发射时刻
	GPSweek t_sig = eph->getShootime(t, psr);
	//用信号发射时刻算卫星位置速度
	result = eph->getSatXvt(t_sig, xvt);
	//地球自转改正
	xvt.pos = eph->getEarthRotCorr((t - t_sig).toSec(), xvt.pos);
	result = 1;
	return result;
}
int EphStore::getXvt(const SatID &sat, GPSweek t, PVT &xvt)const {
	int result = 0;
	const shared_ptr<Ephemeris> eph = findEph(sat, t);
	if (!eph)
		return result;
	result = eph->getSatXvt(t, xvt);
	return result;
}
void EphStore::addEphemeris(const Ephemeris& eph) {
	//if this sat did not exist in the map, create one
	if (ephMap.find(eph.satid) == ephMap.end()) {
		TimeEphMap newTimemap;
		ephMap[eph.satid] = newTimemap;
	}
	TimeEphMap& tem = ephMap[eph.satid];
	////if there is no time map in the sat's map
	//if (tem.size() == 0) {
	//	tem.insert(make_pair(eph.toe, eph.clone()));
	//	return;
	//}
	//if this toe has already existed
	TimeEphMap::iterator item = tem.find(eph.toe);
	if (item == tem.end()) {
		tem.insert(make_pair(eph.toe, eph.clone()));
		return;
	}
	else {
		item->second = eph.clone();
		return;
	}
	tem.insert(make_pair(eph.toe, eph.clone()));
}
const shared_ptr<Ephemeris> EphStore::findEph(const SatID& sat, const GPSweek& t)const {
	if (ephMap.find(sat) == ephMap.end())
		return NULL;
	TimeEphMap tem;
	getTimeEphMap(sat, tem);
	TimeEphMap::const_iterator cit = tem.lower_bound(t);
	if (cit == tem.cbegin()) {
		if (cit->second->isValid(t))
			return cit->second;
		return NULL;
	}
	if (cit == tem.cend()) {
		TimeEphMap::const_reverse_iterator crit = tem.crbegin();
		if (crit->second->isValid(t))
			return crit->second;
		return NULL;
	}
	if (cit->second->isValid(t)) {
		return cit->second;
	}
	cit--;
	if (cit->second->isValid(t)) {
		return cit->second;
	}
	return NULL;
}
bool EphStore::getTimeEphMap(const SatID& sat, EphStore::TimeEphMap& tem)const {
	EphMap::const_iterator cit = ephMap.find(sat);
	if (cit == ephMap.cend()) {
		return false;
	}
	tem = cit->second;
	return true;
}
