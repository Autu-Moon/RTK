#include"PVT.h"
#include<string>
#include<initializer_list>
#pragma once
const double M_PI = 3.14159265358979;				//圆周率


double Deg2Rad(const double& deg);						//角度转弧度
double Rad2Deg(const double& rad);						//弧度转角度
double DCM2Deg(const Triple& dcm);                      //度分秒转角度
double DCM2Rad(const Triple& dcm);                      //度分秒转弧度

double Norm2(std::initializer_list<double> num);

bool isspace(const std::string& str);
template<typename T>
T str2num(const std::string& str){	//string to num
	if (isspace(str))				//if str is space
		return 0;
	size_t D = str.find("D");
	if (D == string::npos) {
		return stod(str);
	}
	else {
		return stod(str.substr(0, D))*
			pow(10, stoi(str.substr(D + 1, str.size() - (D + 1))));
	}
}					

template<typename T> struct map_init_helper
{
	T& data;
	map_init_helper(T& d) : data(d) {}
	map_init_helper& operator() (typename T::key_type const& key, typename T::mapped_type const& value)
	{
		data[key] = value;
		return *this;
	}
};
template<typename T> map_init_helper<T> map_init(T& item)
{
	return map_init_helper<T>(item);
}
