#include"TimeTrans.h"
namespace TimeTrans {
	//标准时到儒略日
	void CommonToMJD(const CommonTime &common, MJuliaday &mjd)
	{
		short year = 0, month = 0;
		if (common.month <= 2) {
			year = common.year - 1;
			month = common.month + 12;
		}
		else {
			year = common.year;
			month = common.month;
		}
		//整数与小数部分分开计算
		mjd.Days = floor(365.25*year) + floor(30.6001*(month + 1)) + common.day + (1720981.5 - 2400000.5);
		mjd.Fracdays = (common.hour + common.minute / 60.0 + common.second / 3600.0) / 24;
	}
	//儒略日到标准时
	void MJDToCommonTime(const MJuliaday &mjd, CommonTime &common)
	{
		//2400000.5即mjd到jd
		unsigned int a = floor(mjd.Days + mjd.Fracdays + 0.5 + 2400000.5);
		unsigned int b = a + 1537;
		unsigned int c = floor((b - 122.1) / 365.25);
		unsigned int d = floor(365.25*c);
		unsigned int e = floor((b - d) / 30.6001);
		double day = b - d - floor(30.6001*e) + fmod(mjd.Fracdays + 0.5 + 2400000.5, 1);
		common.month = e - 1 - 12 * floor(e / 14.0);
		common.year = c - 4715 - floor((7 + common.month) / 10.0);
		common.day = floor(day);
		common.hour = floor(fmod(day, 1) * 24);
		common.minute = floor(fmod(day * 24, 1) * 60);
		common.second = fmod(day * 24 * 60, 1) * 60;
	}

	//commontime to GPSweek
	void CommonToGPSweek(const CommonTime& common, GPSweek& gps) {
		MJuliaday mjd;
		CommonToMJD(common, mjd);
		MJDToGPSweek(mjd, gps);
	}
	//儒略日到GPS周
	void MJDToGPSweek(const MJuliaday &mjd, GPSweek &gps)
	{
		gps.week = floor((mjd.Days - 44244) / 7.0);
		gps.weeksec = (mjd.Days + mjd.Fracdays - 44244 - gps.week * 7) * 86400;
	}
	//GPS周到儒略日
	void GPSweekToMJD(const GPSweek &gps, MJuliaday &mjd)
	{
		mjd.Days = 44244 + gps.week * 7 + floor(gps.weeksec / 86400);
		mjd.Fracdays = fmod(gps.weeksec, 86400) / 86400.0;
	}
	//GPSweek to CommonTime
	void GPSweekToCommon(const GPSweek &gps,CommonTime& ct) {
		MJuliaday mjd;
		GPSweekToMJD(gps, mjd);
		MJDToCommonTime(mjd, ct);
	}


	//commontime to BDSweek
	void CommonToBDSweek(const CommonTime& common, GPSweek& bds) {
		MJuliaday mjd;
		CommonToMJD(common, mjd);
		MJDToBDSweek(mjd, bds);
	}
	//儒略日到BDS周
	void MJDToBDSweek(const MJuliaday &mjd, GPSweek &bds)
	{
		bds.week = floor((mjd.Days - 53736) / 7.0);
		bds.weeksec = (mjd.Days + mjd.Fracdays - 53736 - bds.week * 7) * 86400;
	}
	//BDS周到儒略日
	void BDSweekToMJD(const GPSweek &bds, MJuliaday &mjd)
	{
		mjd.Days = 53736 + bds.week * 7 + floor(bds.weeksec / 86400);
		mjd.Fracdays = fmod(bds.weeksec, 86400) / 86400.0;
	}
	//BDSweek to CommonTime
	void BDSweekToCommon(const GPSweek &bds, CommonTime& ct) {
		MJuliaday mjd;
		BDSweekToMJD(bds, mjd);
		MJDToCommonTime(mjd, ct);
	}


	void BDST2GPST(const GPSweek &bds, GPSweek &gps)
	{
		gps = bds + GPSweek(1356, 14);
	}
	void GPST2BDST(const GPSweek &gps, GPSweek &bds) {
		GPSweek gap(1356, 14);
		bds = gps - gap;
	}
	int getLeapsec(const MJuliaday& mjd)
	{
		int leapsec = 10;
		double JD = mjd.Days;
		if (JD > 41498) leapsec++;
		if (JD > 41682) leapsec++;
		if (JD > 42047) leapsec++;
		if (JD > 42412) leapsec++;
		if (JD > 42777) leapsec++;
		if (JD > 43143) leapsec++;
		if (JD > 43508) leapsec++;
		if (JD > 43873) leapsec++;
		if (JD > 44238) leapsec++;
		if (JD > 44785) leapsec++;
		if (JD > 45150) leapsec++;
		if (JD > 45515) leapsec++;
		if (JD > 46246) leapsec++;
		if (JD > 47160) leapsec++;
		if (JD > 47891) leapsec++;
		if (JD > 48256) leapsec++;
		if (JD > 48803) leapsec++;
		if (JD > 49168) leapsec++;
		if (JD > 49533) leapsec++;
		if (JD > 50082) leapsec++;
		if (JD > 50629) leapsec++;
		if (JD > 51178) leapsec++;
		if (JD > 53735) leapsec++;
		if (JD > 54831) leapsec++;
		if (JD > 56108) leapsec++;
		if (JD > 57203) leapsec++;
		if (JD > 57753) leapsec++;
		return leapsec;
	}
	void UTC2GPST(const CommonTime &common, GPSweek& gps) {
		MJuliaday mjd_utc;
		CommonToMJD(common, mjd_utc);
		CommonToGPSweek(common, gps);
		int leapsec = getLeapsec(mjd_utc);
		gps = gps + leapsec + TAI2GPST;
	}
}