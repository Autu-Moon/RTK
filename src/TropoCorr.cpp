#include"TropoCorr.h"
double TropoCorr::Hopefield(double elev, double H)
{
	double RH = hofElem.RH0*exp(-6.396e-4*(H - hofElem.H0));
	double p = hofElem.P0 * pow((1 - 2.26e-5*(H - hofElem.H0)), 5.225);
	double T = hofElem.T0 - 0.0065*(H - hofElem.H0);
	double e = RH * exp(-37.2465 + 0.213166*T - 2.56908e-4*T*T);
	double h_w = 11000.0;
	double h_d = 40136.0 + 148.72*(hofElem.T0 - 273.16);
	double K_w = 155.2e-7 * 4810 * e * (h_w - H) / (T*T);
	double K_d = 155.2e-7 * p * (h_d - H) / T;
	elev = elev * 180 / M_PI;									//高度角弧度转角度
	double trocorr = K_d / sin(sqrt(elev*elev + 6.25)*M_PI / 180)
		+ K_w / sin(sqrt(elev*elev + 2.25)*M_PI / 180);
	return trocorr;
}
