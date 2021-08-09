#pragma once
#include<iostream>
#include<math.h>
#include"MathUtils.h"

struct Hopefield_Element
{
	//霍普菲尔德模型标准气象元素
	double H0;						//海平面 m
	double T0;						//温度   热力学温度
	double P0;						//气压   mbar
	double RH0;						//相对湿度	
	//默认参数
	Hopefield_Element() {
		H0 = 0;
		T0 = 15 + 273.16;
		P0 = 1013.25;
		RH0 = 0.5;
	}
};
class TropoCorr
{
public:
	double Hopefield(double elev, double H);
private:
	Hopefield_Element hofElem;
};
