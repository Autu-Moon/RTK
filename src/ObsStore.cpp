#include"ObsStore.h"
void ObsStore::addObs(const SatID& sat, const ObsData& obsdata) {
	EpochObsMap::iterator it = EpobsMap.find(sat);
	if (it == EpobsMap.end())
		EpobsMap.insert(make_pair(sat, obsdata));
	else
		it->second = obsdata;
}
void ObsStore::clear() {
	EpobsMap.clear();
}
void ObsStore::dump(std::ostream& os)const {
	os << std::right;
	os.setf(ios::fixed, ios::floatfield);
	os << epoch << endl;
	for (auto sat : EpobsMap) {
		os << sat.first;
		os.precision(3);
		for (size_t i = 0; i < OBSOCCUPIED; i++) {
			os << std::setw(15) << sat.second.Psr[i];
			os << std::setw(15) << sat.second.Adr[i];
			os << std::setw(15) << sat.second.Dop[i];
			os << std::setw(15) << sat.second.CNR[i];
		}
		os << endl;
	}
}