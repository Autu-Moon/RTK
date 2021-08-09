#pragma once
#include<math.h>
#include<iostream>
struct CommonTime	//通用时
{
	unsigned short year;
	unsigned short month;
	unsigned short day;
	unsigned short hour;
	unsigned short minute;
	double second;
	CommonTime(unsigned short year, unsigned short month, unsigned short day, 
		unsigned short hour, unsigned short minute, double second):
		year(year), month(month),day(day),hour(hour),minute(minute),second(second){	}
	CommonTime() = default;
};
struct MJuliaday	//简化儒略日
{
	//整数部分与小数部分分开储存
	unsigned int Days;//范围0~65535
	double Fracdays;
};
class GPSweek		//周+周秒
{
public:
	int week;
	double weeksec;
	GPSweek() {};
	GPSweek(unsigned short w, double s) :week(w), weeksec(s) {};
	GPSweek(const GPSweek &g) {
		this->week = g.week;
		this->weeksec = g.weeksec;
	}
	GPSweek& operator=(const GPSweek& right)
	{
		this->week = right.week;
		this->weeksec = right.weeksec;
		return *this;
	}
	bool operator==(const GPSweek& right)const {
		return((*this - right).toSec() < 1e-7 &&
			(*this - right).toSec() > -1e-7);
	}
	bool operator!=(const GPSweek& right)const {
		return(!(*this==right));
	}
	bool operator<(const GPSweek& right)const
	{
		if (this->week == right.week)
			return(this->weeksec < right.weeksec);
		return(this->week < right.week);
	}
	bool operator>(const GPSweek& right)const
	{
		if (this->week == right.week)
			return((this->weeksec > right.weeksec));
		return(this->week > right.week);
	}
	bool operator<=(const GPSweek& right)const
	{
		if (this->week == right.week)
			return(this->weeksec <= right.weeksec);
		return(this->week <= right.week);
	}
	/*double operator-(const GPSweek& right)const
	{
		return ((this->week - right.week) * 86400 * 7 +
			(this->weeksec - right.weeksec));
	}*/
	GPSweek operator-(double right)const {
		GPSweek temp;
		if (this->weeksec >= right) {
			temp.week = this->week;
			temp.weeksec = this->weeksec - right;

		}
		else {
			temp.week = this->week - 1;
			temp.weeksec = 604800 + this->weeksec - right;
		}
		return temp;
	}
	GPSweek operator+(double right) const {
		GPSweek t;
		if (this->weeksec + right < 604800) {
			t.week = this->week;
			t.weeksec = this->weeksec + right;

		}
		else {
			t.week = this->week + 1;
			t.weeksec = this->weeksec + right - 604800;
		}
		return t;
	}
	GPSweek operator+(const GPSweek& right) const {
		GPSweek t(this->week + right.week, this->weeksec + right.weeksec);
		if (t.weeksec >= 604800) {
			t.weeksec = t.weeksec - 604800;
			t.week++;
		}
		return t;
	}
	GPSweek operator-(const GPSweek& right) const {
		GPSweek t(this->week - right.week, this->weeksec - right.weeksec);
		if (t.weeksec < 0) {
			t.weeksec = t.weeksec + 604800;
			t.week--;
		}
		return t;
	}
	friend 	std::ostream& operator<<(std::ostream& out, const GPSweek& gps);
	double toSec() {
		return(this->week * 604800 + this->weeksec);
	}
};

inline std::ostream& operator<<(std::ostream& out, const GPSweek& gps)
{
	out << gps.week << "  " << gps.weeksec;
	return out;
}
namespace TimeTrans {
	void CommonToMJD(const CommonTime &common, MJuliaday &mjd);
	void MJDToCommonTime(const MJuliaday &mjd, CommonTime &common);
	/*
	* @Description: TimeTrans in GPST
	*/
	void MJDToGPSweek(const MJuliaday &mjd, GPSweek &gps);
	void CommonToGPSweek(const CommonTime& common, GPSweek& gps);
	void GPSweekToMJD(const GPSweek &gps, MJuliaday &mjd);
	void GPSweekToCommon(const GPSweek &gps, CommonTime &ct);
	/*
	* @Description: TimeTrans in BDT
	*/
	void MJDToBDSweek(const MJuliaday &mjd, GPSweek &gps);
	void CommonToBDSweek(const CommonTime& common, GPSweek& gps);
	void BDSweekToMJD(const GPSweek &gps, MJuliaday &mjd);
	void BDSweekToCommon(const GPSweek &gps, CommonTime &ct);
	/*
	* @Description: TimeTrans between GPST and BDT
	*/
	void BDST2GPST(const GPSweek &bds, GPSweek &gps);
	void GPST2BDST(const GPSweek &bds, GPSweek &gps);


	/*
	* @Description: give leapsecond in UTC by TAI
	* @Param: mjd of UTC
	* @return: leapsecond
	*/
	int getLeapsec(const MJuliaday &mjd);
	const int TAI2GPST = -19;			//TAI is always ahead of GPST by 19s
	/*
	* @Description: trans UTC to GPST
	* @Param: commontime of UTC
	* @return: gpsweek of gpst
	*/
	void UTC2GPST(const CommonTime &common, GPSweek& gps);
}

