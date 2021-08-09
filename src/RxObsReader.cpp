#include "RxObsReader.h"
#include<stdlib.h>

bool Rx304ObsReader::readHead(std::ifstream& fin, ObsStore& obsStore) {
	std::string line;
	std::getline(fin, line);
	/*if (line.substr(0, 9).find("3.04")== std::string::npos) {
		std::cerr << "Error in Rinex Edition, expecting 3.04" << std::endl;
		return false;
	}*/
	while (!fin.eof())
	{
		//判断
		std::getline(fin, line);
		vector<ObsCode>* typevec(NULL);
		if (line.find(TYPE, 0) != std::string::npos)
		{
			//store Observation type data
			//系统
			std::string sys = line.substr(0, 1);
			if (sys == "G") {
				typevec = &CodeVec[GPS];
			}
			else if (sys == "C") {
				typevec = &CodeVec[BDS];
			}
			else
				continue;
			//观测类型数
			int typeNumber = stoi(line.substr(4, 2));
			//根据观测类型数储存观测类型名
			if (typeNumber > 13)
			{
				for (int i = 0; i < 13; i++)
				{
					ObsCode type = str2ObsCode(line.substr(7 + 4 * i, 3));
					if (type) {
						typevec->push_back(type);
					}
				}
				std::getline(fin, line);
				for (int i = 0; i < typeNumber - 13; i++)
				{
					ObsCode type = str2ObsCode(line.substr(7 + 4 * i, 3));
					if (type) {
						typevec->push_back(type);
					}
				}
			}
			else if (typeNumber > 26)
			{
				for (int i = 0; i < 13; i++)
				{
					ObsCode type = str2ObsCode(line.substr(7 + 4 * i, 3));
					if (type) {
						typevec->push_back(type);
					}
				}
				std::getline(fin, line);
				for (int i = 0; i < typeNumber - 13; i++)
				{
					ObsCode type = str2ObsCode(line.substr(7 + 4 * i, 3));
					if (type) {
						typevec->push_back(type);
					}
				}
				std::getline(fin, line);
				for (int i = 0; i < typeNumber - 26; i++)
				{
					ObsCode type = str2ObsCode(line.substr(7 + 4 * i, 3));
					if (type) {
						typevec->push_back(type);
					}
				}
			}
			else if (typeNumber <= 13)
			{
				for (int i = 0; i < 13; i++)
				{
					ObsCode type = str2ObsCode(line.substr(7 + 4 * i, 3));
					if (type) {
						typevec->push_back(type);
					}
				}
			}
		}
		if (line.find(Interval, 0) != line.npos)
		{
			obsStore.obsid.interval = stoi(line.substr(0, 20));
		}
		if (line.find(EOH, 0) != line.npos)
			break;
	}
	return true;
}

int  Rx304ObsReader::readData(std::ifstream& fin, ObsStore& obsStore) {
	obsStore.clear();
	epobsMap.clear();
	string line;
	while (!fin.eof())
	{
		std::getline(fin, line);
		//find the start of epoch
		if (line.substr(0, 1) != EpochMark)
			continue;
		//epoch info
		unsigned short year = (unsigned short)stoi(line.substr(2, 4));
		unsigned short month = (unsigned short)stoi(line.substr(7, 2));
		unsigned short day = (unsigned short)stoi(line.substr(10, 2));
		unsigned short hour = (unsigned short)stoi(line.substr(13, 2));
		unsigned short minute = (unsigned short)stoi(line.substr(16, 2));
		double second = stod(line.substr(18, 11));
		TimeTrans::CommonToGPSweek(CommonTime(year, month, day, hour, minute, second), obsStore.epoch);
		//num of satellite
		int satNum = stoi(line.substr(32, 3));
		//obs data
		for (int i = 0; i < satNum; i++)
		{
			SYSTEM mark;
			std::getline(fin, line);
			if (line.substr(0, 1) == "G") {					//GPS
				mark = GPS;
			}
			else if (line.substr(0, 1) == "C") {			//BDS
				mark = BDS;
			}
			else
				continue;
			fillblank(line);
			SatID  sat(mark, stoi(line.substr(1, 2)));
			for (size_t j = 0; j < CodeVec[mark].size(); j++)
			{
				Observation obs;
				obs.obs = str2num<double>(line.substr(3 + 16 * j, 14));
				obs.LLI = str2num<short>(line.substr(17 + 16 * j, 1));
				obs.SSI = str2num<short>(line.substr(18 + 16 * j, 1));
				if (obs.obs == 0 && obs.LLI == 0 && obs.SSI == 0)
					obs.valid = 0;
				else
					obs.valid = 1;
				addObservation(sat, CodeVec[mark][j], obs);
			}

		}
		toStored(obsStore);
		break;//only for a single epoch
	}
	return 1;
}

void Rx304ObsReader::addObservation(const SatID& sat, const ObsCode& code, const Observation& obs) {
	if (epobsMap.find(sat) == epobsMap.end()) {	//if no sat
		ObsCodeMap newCodeMap;
		epobsMap.insert(make_pair(sat, newCodeMap));
	}
	ObsCodeMap& tim = epobsMap[sat];
	ObsCodeMap::iterator itim = tim.find(code);
	if (itim == tim.end()) {			//if no typeID
		tim.insert(make_pair(code, obs));
		return;
	}
	itim->second = obs;
}

void Rx304ObsReader::invalid() {
	for (auto& sat : epobsMap) {
		for (auto& obs : sat.second) {
			obs.second.valid = 0;
		}
	}
}

void Rx304ObsReader::toStored(ObsStore& obsStore) {
	bool _used = false;
	ObsCodeMap::iterator itim;
	ObsData tobs;
	obsStore.clear();
	for (auto sat : epobsMap) {
		ObsCodeMap& typemap = sat.second;
		tobs.reset();
		_used = false;
		for (auto code : expCodeVec[sat.first.sys]) {
			//make sure this obs exist
			itim = typemap.find(code);
			if (itim == typemap.end())
				continue;
			if (itim->second.valid == 0)
				continue;
			//location of this obs in ObsData based on freq
			size_t fband = getLocation(code, sat.first.sys);
			//psr adr cnr dop
			if (obscodeString[code].substr(0, 1) == obscodeString[_C]) {
				tobs.Psr[fband] = itim->second.obs;
				_used = true;
			}
			if (obscodeString[code].substr(0, 1) == obscodeString[_L]) {
				tobs.Adr[fband] = itim->second.obs;
				_used = true;
			}
			if (obscodeString[code].substr(0, 1) == obscodeString[_D]) {
				tobs.Dop[fband] = itim->second.obs;
				_used = true;
			}
			if (obscodeString[code].substr(0, 1) == obscodeString[_S]) {
				tobs.CNR[fband] = itim->second.obs;
				_used = true;
			}
		}
		if (_used == true) {
			obsStore.addObs(sat.first, tobs);
		}
	}

}

void Rx304ObsReader::setExpectCode(SYSTEM sys,initializer_list<string> code) {
	for (auto _c : code) {
		expCodeVec[sys].push_back(str2ObsCode(_c));
	}
}

void Rx304ObsReader::fillblank(string& str) {
	str = str + "                                                                                                                                                                                                    ";
}