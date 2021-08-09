#include "UserObsStore.h"

void UserObsStore::setStore(ObsStore& obsStore,const EphStore& ephStore) {
	clear();
	PVT pvt;
	UserObs data;
	shared_ptr<Ephemeris> eph(NULL);
	this->epoch = obsStore.epoch;
	//观测值逐个迭代
	for(auto itobs:obsStore.EpobsMap) {
		if (!ephStore.getXvt(itobs.first, itobs.second.Psr[0], obsStore.epoch, data.pvt)) {
			continue;
		}
		eph = ephStore.findEph(itobs.first, obsStore.epoch);
		if (!eph) {
			continue;
		}
		memcpy(data.Tgd, eph->Tgd, EPHTGDLENGTH * sizeof(double));
		memcpy(data.Psr, itobs.second.Psr, OBSERVATIONLENGTH * sizeof(double));
		memcpy(data.Adr, itobs.second.Adr, OBSERVATIONLENGTH * sizeof(double));
		memcpy(data.Dop, itobs.second.Dop, OBSERVATIONLENGTH * sizeof(double));
		memcpy(data.CNR, itobs.second.CNR, OBSERVATIONLENGTH * sizeof(double));
		addUserObs(itobs.first, data);
	}
	//obsStore.clear();
}
void UserObsStore::addUserObs(const SatID& sat, const UserObs& data) {
	UserObsMap::iterator isom = userobsMap.find(sat);
	if (isom == userobsMap.end()) {
		userobsMap.insert(make_pair(sat, data));
	}
	else
		isom->second = data;
}
void UserObsStore::dump(std::ostream &os)const {
	os << std::right;
	os.setf(ios::fixed, ios::floatfield);
	os << epoch << "   " << userobsMap.size() << endl;
	for (auto sat : userobsMap) {
		os << sat.first;
		os.precision(4);
		os << std::setw(15) << sat.second.pvt.pos;
		//os << std::setw(13) << sat.second.pvt.vel;
		os.precision(10);
		os << std::setw(13) << sat.second.pvt.clockbias;
		//os << std::setw(13) << sat.second.pvt.clockdrift;
		os.precision(4);
		for (size_t i = 0; i < OBSOCCUPIED; i++) {
			os << std::setw(15) << sat.second.Psr[i];
			os << std::setw(15) << sat.second.Adr[i];

		}
		os << std::setw(15) << sat.second.IFComb;
		os << endl;
	}
}
//bool UserObsStore::setComb(CombineObs Comb) {
//	if (Comb == IonoFree) {
//		setIFComb();
//	}
//	return true;
//}
void UserObsStore::setIFComb(initializer_list<FreqBand> _band) {
	vector<double> freq;
	vector<FreqBand> band;
	for (auto _b : _band) {
		freq.push_back(getFreq(_b));
		band.push_back(_b);
	}
	double ifobs;
	if (band[0] < 3) {//GPS
		//L1 L2
		for (UserObsMap::iterator it = userobsMap.begin();
			it != userobsMap.end(); it++)										//逐卫星迭代
		{
			if (it->first.sys != GPS)continue;
			//计算消电离层组合观测值
			if (it->second.Psr[band[0]] != OBSEMPTY && it->second.Psr[band[1]] != OBSEMPTY) {
				ifobs = freq[0] * freq[0] / (freq[0] * freq[0] - freq[1] * freq[1])*it->second.Psr[band[0]] -
					freq[1] * freq[1] / (freq[0] * freq[0] - freq[1] * freq[1])*it->second.Psr[band[1]];
				it->second.IFComb = ifobs;
			}
			else
				continue;
		}
	}
	else if (band[0] >= 3) {//BDS
		//B1I B3I
		for (UserObsMap::iterator it = userobsMap.begin();
			it != userobsMap.end(); it++) {								//逐卫星迭代
			//计算消电离层组合观测值		
			if (it->first.sys != BDS)
				continue;
			if (it->second.Psr[0] != OBSEMPTY && it->second.Psr[1] != OBSEMPTY) {
				ifobs = freq[0] * freq[0] / (freq[0] * freq[0] - freq[1] * freq[1])*it->second.Psr[0] -
					freq[1] * freq[1] / (freq[0] * freq[0] - freq[1] * freq[1])*it->second.Psr[1];
				//电离层的群延
				ifobs -= V_light * freq[0] * freq[0] * it->second.Tgd[0] / (freq[0] * freq[0] - freq[1] * freq[1]);
				it->second.IFComb = ifobs;
			}
			else
				continue;
		}
	}
}
void UserObsStore::setElev(const Triple& pos) {
	CoorTrans _t;
	for (auto &iubs : userobsMap) {
		if (!_t.XYZToNEU(iubs.second.pvt.pos, pos, iubs.second.elev)) {
			continue;
		}
	}
}
void UserObsStore::setCutOffElev(double cutoff) {
	this->CutOffElev = Deg2Rad(cutoff);
	for (UserObsMap::iterator it = userobsMap.begin(); it != userobsMap.end(); ) {
		if (it->second.elev < CutOffElev) {
			userobsMap.erase(it++);
			continue;
		}
		it++;
	}
}
void UserObsStore::setBanSat(initializer_list<SatID> banSat) {
	for (auto sat : banSat) {
		userobsMap.erase(sat);
	}
}
void UserObsStore::setBanBDSGEO() {
	setBanSat({ SatID(BDS,1),SatID(BDS,2),SatID(BDS,3),SatID(BDS,4),SatID(BDS,5),
		SatID(BDS,59),SatID(BDS,60),SatID(BDS,61) });
}
void UserObsStore::setBanBDSIGSO() {
	setBanSat({ SatID(BDS,6),SatID(BDS,7),SatID(BDS,8),SatID(BDS,9),SatID(BDS,10),
		SatID(BDS,13),SatID(BDS,16),SatID(BDS,38),SatID(BDS,39),SatID(BDS,40) });
}
void UserObsStore::setBanBDS2() {
	//GEO
	setBanSat({ SatID(BDS,1),SatID(BDS,2),SatID(BDS,3),SatID(BDS,4),SatID(BDS,5) });
	//IGSO
	setBanSat({ SatID(BDS,6),SatID(BDS,7),SatID(BDS,8),SatID(BDS,9),SatID(BDS,10),
			SatID(BDS,13),SatID(BDS,16) });
	//MEO
	setBanSat({ SatID(BDS,11),SatID(BDS,12),SatID{BDS,14} });
}
void UserObsStore::setBanBDS3() {
	//IGSO
	setBanSat({ SatID(BDS,38),SatID(BDS,39),SatID(BDS,40) });
	//MEO
	setBanSat({ SatID(BDS,19), SatID(BDS, 20), SatID(BDS, 21),SatID(BDS,22),SatID(BDS,23),SatID(BDS,24),
		SatID(BDS,25),SatID(BDS,26),SatID(BDS,27),SatID(BDS,28),SatID(BDS,29),SatID(BDS,30),SatID(BDS,32),
		SatID(BDS,33),SatID(BDS,34),SatID(BDS,35),SatID(BDS,36),SatID(BDS,37),SatID(BDS,41),SatID(BDS,42),
		SatID(BDS,43),SatID(BDS,44),SatID(BDS,45),SatID(BDS,46) });
	//GEO
	setBanSat({ SatID(BDS,59),SatID(BDS,60),SatID(BDS,61) });
}
void UserObsStore::setBanBDS() {
	for (UserObsMap::iterator it = userobsMap.begin(); it != userobsMap.end();) {
		if (it->first.sys == BDS)
			userobsMap.erase(it++);
		else
			it++;
	}
}
void UserObsStore::setBanGPS() {
	for (UserObsMap::iterator it = userobsMap.begin(); it != userobsMap.end();) {
		if (it->first.sys == GPS)
			userobsMap.erase(it++);
		else
			it++;
	}
}
void UserObsStore::setFreq(initializer_list<FreqBand> band) {
	vector<size_t> SYSBand[SYScount];
	bool valid = false;
	for (auto _b : band) {
		if (_b >= L1 && _b <= L5) {
			SYSBand[GPS].push_back(getLocation(_b));
		}
		else if (_b >= B1I && _b <= B3I) {
			SYSBand[BDS].push_back(getLocation(_b));
		}
	}
	for (UserObsMap::iterator itobs = userobsMap.begin(); itobs != userobsMap.end();) {
		valid = true;
		if (itobs->first.sys > BDS) {
			valid = false;
		}
		for (auto _b : SYSBand[itobs->first.sys]) {
			if (itobs->second.Psr[_b] == OBSEMPTY) {
				valid = false;
				break;
			}
			if (itobs->second.Adr[_b] == OBSEMPTY) {
				valid = false;
				break;
			}
		}
		if (valid == false) {
			userobsMap.erase(itobs++);
		}
		else
			itobs++;
	}
}

SatID UserObsStore::getBest(SYSTEM sys)const {
	SatID sat;
	double elev = 0;
	for (auto iubs : userobsMap) {
		if (iubs.first.sys != sys)
			continue;
		//make sure no lack of obs
		if (iubs.second.elev > elev && !iubs.second.isLack() && iubs.second.best == 1) {
			elev = iubs.second.elev;
			sat = iubs.first;
		}
		else
			continue;
	}
	return sat;
}
void UserObsStore::clear() {
	userobsMap.clear();
}
int UserObsStore::size()const {
	return this->userobsMap.size();
}
