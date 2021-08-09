#include"RefFrame.h"
RefFrame::RefFrame(RefFrameType r) {
	ref = r;
	switch (ref)
	{
	case WGS84:
		a = 6.3781370e6;
		f = 1 / 298.257223563;
		GM = 3.986004418e14;
		g_a = 9.7803267714;
		U_0 = 62636860.8497;
		esq = 0.00669437999013;
		// esq = 2 * f - f * f;
		break;
	case CGCS2000:
		a = 6378137;
		f = 1 / 298.257222101;
		GM = 3.986004418e14;
		esq = 2 * f - f * f;
		break;
	case GTRF:
		a = 6.37813655e6;
		f = 1 / 298.25769;
		GM = 3.986004415e14;
		esq = 2 * f - f * f;
		break;
	case ITRF:
		a = 6.3781366e6;
		f = 1 / 298.25642;
		GM = 3.986004418e14;
		U_0 = 6.26368560e7;
		g_a = 9.7803278;
		esq = 2 * f - f * f;
		break;
	default:
		break;
	}
	b = (1 - f) * a;
}
void RefFrame::setFrame(RefFrameType r) {
	ref = r;
	switch (ref)
	{
	case WGS84:
		a = 6.3781370e6;
		f = 1 / 298.257223563;
		GM = 3.986004418e14;
		g_a = 9.7803267714;
		U_0 = 62636860.8497;
		break;
	case CGCS2000:
		a = 6378137;
		f = 1 / 298.257222101;
		GM = 3.986004418e14;
	case GTRF:
		a = 6.37813655e6;
		f = 1 / 298.25769;
		GM = 3.986004415e14;
		break;
	case ITRF:
		a = 6.3781366e6;
		f = 1 / 298.25642;
		GM = 3.986004418e14;
		U_0 = 6.26368560e7;
		g_a = 9.7803278;
		break;
	default:
		break;
	}
	esq = 2 * f - f * f;
	b = (1 - f) * a;
}
double RefFrame::computeR_N(double B) {
	return (a / sqrt(1. - esq * sin(B) * sin(B)));
}
double RefFrame::computeR_M(double B) {
	return (a * (1 - esq)) /
		pow((1. - esq * sin(B) * sin(B)), 1.5);
}