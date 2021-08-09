#include"MathUtils.h"
double Deg2Rad(const double& deg)						//角度转弧度
{
	return(deg / 180.0*M_PI);
}
double Rad2Deg(const double& rad)						//弧度转角度
{
	return(rad / M_PI * 180.0);
}
double DCM2Deg(const Triple& dcm)
{
	return(dcm[0] + dcm[1] / 60.0 + dcm[2] / 3600.0);
}
double DCM2Rad(const Triple& dcm)
{
	return(Deg2Rad(DCM2Deg(dcm)));
}
double Norm2(std::initializer_list<double> num) {
	double res = 0;
	for (auto _n : num) {
		res += _n * _n;
	}
	res = sqrt(res);
	return res;
}

bool isspace(const std::string& str)
{
	for (auto s : str) {
		if (s == '\t' || s == '\n' || s == ' ')
			continue;
		return false;
	}
	return true;
}