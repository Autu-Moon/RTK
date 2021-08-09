#pragma once
#include<valarray>
#include<iostream>
class Triple
{
public:
	std::valarray<double> thearray;
public:
	Triple() :thearray(3) {};
	Triple(const Triple& right) :thearray(right.thearray) {};
	Triple(double a, double b, double c) :thearray(3)
	{
		thearray[0] = a;
		thearray[1] = b;
		thearray[2] = c;
	}
	Triple operator-(const Triple& right)const {
		Triple t;
		t[0] = this->thearray[0] - right[0];
		t[1] = this->thearray[1] - right[1];
		t[2] = this->thearray[2] - right[2];
		return t;
	}
	Triple operator+(const Triple& right)const {
		Triple t;
		t[0] = this->thearray[0] + right[0];
		t[1] = this->thearray[1] + right[1];
		t[2] = this->thearray[2] + right[2];
		return t;
	}
	double operator[](const size_t index) const
	{
		return thearray[index];
	}
	double& operator[](const size_t index)
	{
		return thearray[index];
	}
	Triple& operator= (const Triple &right)
	{
		this->thearray[0] = right[0];
		this->thearray[1] = right[1];
		this->thearray[2] = right[2];
		return *this;
	}
	bool operator== (const Triple &right)
	{
		return(thearray[0] == right[0] && thearray[1] == right[1] && thearray[2] == right[2]);
	}
	bool operator!= (const Triple &right)
	{
		return(thearray[0] != right[0] || thearray[1] != right[1] || thearray[2] != right[2]);
	}
	double norm() {
		return sqrt(thearray[0] * thearray[0] +
			thearray[1] * thearray[1] + thearray[2] * thearray[2]);
	}
	friend std::ostream& operator<<(std::ostream &out, const Triple t);
};
inline std::ostream& operator<<(std::ostream &out, const Triple t)
{
	out << "(" << t[0] << " , " << t[1] << " , " << t[2] << ")";
	return out;
}
class PVT
{
public:
	Triple pos;
	Triple vel;
	double clockbias;
	double clockdrift;
	PVT() :pos(0.0, 0.0, 0.0), vel(0.0, 0.0, 0.0) {};
};
