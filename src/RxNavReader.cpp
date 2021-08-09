#include "RxNavReader.h"
int Rx304NavReader::read(std::ifstream& fin, EphStore& ephStore) {
	readHead(fin, ephStore);
	readData(fin, ephStore);
	return 1;
}
int Rx304NavReader::readData(std::ifstream& fin, EphStore& ephStore) {
	string line;	
	GPSEphemeris gpseph;
	BDSEphemeris bdseph;
	GPSweek toc;
	while (!fin.eof()) {
		getline(fin, line);	//1
		if (line.substr(0, 1) == "G") {
			unsigned short year = (unsigned short)stoi(line.substr(4, 4));
			unsigned short month = (unsigned short)stoi(line.substr(8, 3));
			unsigned short day = (unsigned short)stoi(line.substr(11, 3));
			unsigned short hour = (unsigned short)stoi(line.substr(14, 3));
			unsigned short minute = (unsigned short)stoi(line.substr(17, 3));
			double second = stod(line.substr(20, 3));
			TimeTrans::CommonToGPSweek(CommonTime(year, month, day, hour, minute, second), toc);
			gpseph.toc = toc;
			gpseph.satid = SatID(GPS, stoi(line.substr(1, 2)));
			gpseph.af0 = (str2num<double>(line.substr(23, 19)));
			gpseph.af1 = (str2num<double>(line.substr(42, 19)));
			gpseph.af2 = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//2
			gpseph.IODE = (str2num<int>(line.substr(4, 19)));
			gpseph.crs = (str2num<double>(line.substr(23, 19)));
			gpseph.deltaN = (str2num<double>(line.substr(42, 19)));
			gpseph.M0 = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//3
			gpseph.cuc = (str2num<double>(line.substr(4, 19)));
			gpseph.ecc = (str2num<double>(line.substr(23, 19)));
			gpseph.cus = (str2num<double>(line.substr(42, 19)));
			gpseph.rootA = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//4
			gpseph.toe = GPSweek(toc.week, str2num<double>(line.substr(4, 19)));
			gpseph.cic = (str2num<double>(line.substr(23, 19)));
			gpseph.Omega0 = (str2num<double>(line.substr(42, 19)));
			gpseph.cis = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//5
			gpseph.i0 = (str2num<double>(line.substr(4, 19)));
			gpseph.crc = (str2num<double>(line.substr(23, 19)));
			gpseph.omega = (str2num<double>(line.substr(42, 19)));
			gpseph.Omega_dot = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//6
			gpseph.i_dot = (str2num<double>(line.substr(4, 19)));
			getline(fin, line);	//7
			gpseph.Tgd[0] = (str2num<double>(line.substr(42, 19)));
			getline(fin, line);	//8
			gpseph.tow = (str2num<double>(line.substr(4, 19)));
			ephStore.addEphemeris(gpseph);
		}
		else if (line.substr(0, 1) == "C") {
			unsigned short year = (unsigned short)stoi(line.substr(4, 4));
			unsigned short month = (unsigned short)stoi(line.substr(8, 3));
			unsigned short day = (unsigned short)stoi(line.substr(11, 3));
			unsigned short hour = (unsigned short)stoi(line.substr(14, 3));
			unsigned short minute = (unsigned short)stoi(line.substr(17, 3));
			double second = stod(line.substr(20, 3));
			TimeTrans::CommonToBDSweek(CommonTime(year, month, day, hour, minute, second), toc);
			TimeTrans::BDST2GPST(toc, bdseph.toc);				//store as GPST
			bdseph.satid = SatID(BDS, stoi(line.substr(1, 2)));
			bdseph.af0 = (str2num<double>(line.substr(23, 19)));
			bdseph.af1 = (str2num<double>(line.substr(42, 19)));
			bdseph.af2 = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//2
			bdseph.AODE = (str2num<int>(line.substr(4, 19)));
			bdseph.crs = (str2num<double>(line.substr(23, 19)));
			bdseph.deltaN = (str2num<double>(line.substr(42, 19)));
			bdseph.M0 = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//3
			bdseph.cuc = (str2num<double>(line.substr(4, 19)));
			bdseph.ecc = (str2num<double>(line.substr(23, 19)));
			bdseph.cus = (str2num<double>(line.substr(42, 19)));
			bdseph.rootA = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//4
			
			bdseph.toe = GPSweek(toc.week, str2num<double>(line.substr(4, 19)));
			TimeTrans::BDST2GPST(bdseph.toe, bdseph.toe);
			bdseph.cic = (str2num<double>(line.substr(23, 19)));
			bdseph.Omega0 = (str2num<double>(line.substr(42, 19)));
			bdseph.cis = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//5
			bdseph.i0 = (str2num<double>(line.substr(4, 19)));
			bdseph.crc = (str2num<double>(line.substr(23, 19)));
			bdseph.omega = (str2num<double>(line.substr(42, 19)));
			bdseph.Omega_dot = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//6
			bdseph.i_dot = (str2num<double>(line.substr(4, 19)));
			getline(fin, line);	//7
			bdseph.Tgd[0] = (str2num<double>(line.substr(42, 19)));
			bdseph.Tgd[1] = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//8
			bdseph.tow = (str2num<double>(line.substr(4, 19)));
			ephStore.addEphemeris(bdseph);
		}
		else
			continue;
	}
	return 1;

}
bool Rx304NavReader::readHead(std::ifstream& fin, EphStore& ephStore) {
	string line;
	while (!fin.eof())
	{
		getline(fin, line);
		if (line.find("END OF HEADER", 0) != line.npos)
			break;
	}
	return true;
}

int Rx302NavReader::read(std::ifstream& fin, EphStore& ephStore) {
	readHead(fin, ephStore);
	readData(fin, ephStore);
	return 1;
}

int Rx302NavReader::readData(std::ifstream& fin, EphStore& ephStore) {
	string line;
	GPSEphemeris gpseph;
	BDSEphemeris bdseph;
	GPSweek toc;
	while (!fin.eof()) {
		getline(fin, line);	//1
		if (line.substr(0, 1) == "G") {
			unsigned short year = (unsigned short)stoi(line.substr(4, 4));
			unsigned short month = (unsigned short)stoi(line.substr(8, 3));
			unsigned short day = (unsigned short)stoi(line.substr(11, 3));
			unsigned short hour = (unsigned short)stoi(line.substr(14, 3));
			unsigned short minute = (unsigned short)stoi(line.substr(17, 3));
			double second = stod(line.substr(20, 3));
			TimeTrans::CommonToGPSweek(CommonTime(year, month, day, hour, minute, second), toc);
			gpseph.toc = toc;
			gpseph.satid = SatID(GPS, stoi(line.substr(1, 2)));
			gpseph.af0 = (str2num<double>(line.substr(23, 19)));
			gpseph.af1 = (str2num<double>(line.substr(42, 19)));
			gpseph.af2 = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//2
			gpseph.IODE = (str2num<int>(line.substr(4, 19)));
			gpseph.crs = (str2num<double>(line.substr(23, 19)));
			gpseph.deltaN = (str2num<double>(line.substr(42, 19)));
			gpseph.M0 = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//3
			gpseph.cuc = (str2num<double>(line.substr(4, 19)));
			gpseph.ecc = (str2num<double>(line.substr(23, 19)));
			gpseph.cus = (str2num<double>(line.substr(42, 19)));
			gpseph.rootA = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//4
			gpseph.toe = GPSweek(toc.week, str2num<double>(line.substr(4, 19)));
			gpseph.cic = (str2num<double>(line.substr(23, 19)));
			gpseph.Omega0 = (str2num<double>(line.substr(42, 19)));
			gpseph.cis = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//5
			gpseph.i0 = (str2num<double>(line.substr(4, 19)));
			gpseph.crc = (str2num<double>(line.substr(23, 19)));
			gpseph.omega = (str2num<double>(line.substr(42, 19)));
			gpseph.Omega_dot = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//6
			gpseph.i_dot = (str2num<double>(line.substr(4, 19)));
			getline(fin, line);	//7
			gpseph.Tgd[0] = (str2num<double>(line.substr(42, 19)));
			getline(fin, line);	//8
			gpseph.tow = (str2num<double>(line.substr(4, 19)));
			ephStore.addEphemeris(gpseph);
		}
		else if (line.substr(0, 1) == "C") {
			unsigned short year = (unsigned short)stoi(line.substr(4, 4));
			unsigned short month = (unsigned short)stoi(line.substr(8, 3));
			unsigned short day = (unsigned short)stoi(line.substr(11, 3));
			unsigned short hour = (unsigned short)stoi(line.substr(14, 3));
			unsigned short minute = (unsigned short)stoi(line.substr(17, 3));
			double second = str2num<double>(line.substr(20, 3));
			TimeTrans::CommonToBDSweek(CommonTime(year, month, day, hour, minute, second), toc);
			TimeTrans::BDST2GPST(toc, bdseph.toc);				//store as GPST
			bdseph.satid = SatID(BDS, stoi(line.substr(1, 2)));
			bdseph.af0 = (str2num<double>(line.substr(23, 19)));
			bdseph.af1 = (str2num<double>(line.substr(42, 19)));
			bdseph.af2 = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//2
			bdseph.AODE = (str2num<int>(line.substr(4, 19)));
			bdseph.crs = (str2num<double>(line.substr(23, 19)));
			bdseph.deltaN = (str2num<double>(line.substr(42, 19)));
			bdseph.M0 = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//3
			bdseph.cuc = (str2num<double>(line.substr(4, 19)));
			bdseph.ecc = (str2num<double>(line.substr(23, 19)));
			bdseph.cus = (str2num<double>(line.substr(42, 19)));
			bdseph.rootA = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//4

			bdseph.toe = GPSweek(toc.week, str2num<double>(line.substr(4, 19)));
			TimeTrans::BDST2GPST(bdseph.toe, bdseph.toe);
			bdseph.cic = (str2num<double>(line.substr(23, 19)));
			bdseph.Omega0 = (str2num<double>(line.substr(42, 19)));
			bdseph.cis = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//5
			bdseph.i0 = (str2num<double>(line.substr(4, 19)));
			bdseph.crc = (str2num<double>(line.substr(23, 19)));
			bdseph.omega = (str2num<double>(line.substr(42, 19)));
			bdseph.Omega_dot = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//6
			bdseph.i_dot = (str2num<double>(line.substr(4, 19)));
			getline(fin, line);	//7
			bdseph.Tgd[0] = (str2num<double>(line.substr(42, 19)));
			bdseph.Tgd[1] = (str2num<double>(line.substr(61, 19)));
			getline(fin, line);	//8
			bdseph.tow = (str2num<double>(line.substr(4, 19)));
			ephStore.addEphemeris(bdseph);
		}
		else
			continue;
	}
	return 1;

}
bool Rx302NavReader::readHead(std::ifstream& fin, EphStore& ephStore) {
	string line;
	while (!fin.eof())
	{
		getline(fin, line);
		if (line.find("END OF HEADER", 0) != line.npos)
			break;
	}
	return true;
}