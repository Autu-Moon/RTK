#pragma once
#include<iostream>
#include<cmath>

class RefFrame
{
public:
	enum RefFrameType {				//参考框架
		WGS84 = 0,					//GPS
		ITRF,                       //ITRS
		CGCS2000,                   //BDS
		GTRF,                       //Galileo
		PZ90                        //GLONASS
	};
	//-------------field-------------
	RefFrameType    ref;
	double          a;                      //长半轴 unit：m
	double          b;                      //短半轴
	double          f;                      //扁率
	double          esq;					//第一偏心率平方
	double          GM;                     //地球引力常数 unit:m^3/s^2
	double          g_a;                    //赤道正常重力 unit:m/s^2
	double          U_0;                    //椭球正常重力位 unit:m^2/s^2
	//-------------constructor-------------
	RefFrame(RefFrameType r = WGS84);
	RefFrame(const RefFrame& right);
	void setFrame(RefFrameType r);
	//-------------function-------------
	double computeR_N(double B);                  //计算卯酉圈半径N
	double computeR_M(double B);                  //计算子午圈半径M
	//-------------operator-------------
	RefFrame& operator= (const RefFrame &right)
	{
		this->ref = right.ref;
		this->a = right.a;
		this->b = right.b;
		this->f = right.f;
		this->esq = right.esq;
		this->GM = right.GM;
		this->g_a = right.g_a;
		this->U_0 = right.U_0;
		return *this;
	}
};