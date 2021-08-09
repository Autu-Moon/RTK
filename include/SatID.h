#pragma once
#include"SatNavBasic.h"
#include<iostream>
class SatID
{
public:
	SYSTEM sys;
	short svprn;
	SatID(SYSTEM s, short prn)
	{
		sys = s;
		svprn = prn;
	}
	SatID() {
		sys = SYSunknown;
		svprn = -1;
	};
	bool isEmpty()const {
		return((sys == SYSunknown) && (svprn == -1));
	}
	SatID& operator=(const SatID& right)
	{
		this->sys = right.sys;
		this->svprn = right.svprn;
		return *this;
	}
	bool operator==(const SatID& right) const
	{
		return ((sys == right.sys) && (svprn == right.svprn));
	}
	bool operator<(const SatID& right) const
	{
		if (sys == right.sys)
			return (svprn < right.svprn);
		return (sys < right.sys);
	}
	friend std::ostream& operator<< (std::ostream& out, const SatID& id);
	void reset() {
		sys = SYSunknown;
		svprn = -1;
	}
};
inline std::ostream& operator<< (std::ostream& out, const SatID& id)
{
	if (id.sys == GPS) {
		out << "G" << id.svprn;
	}
	if (id.sys == BDS) {
		out << "C" << id.svprn;
	}
	return out;
}
