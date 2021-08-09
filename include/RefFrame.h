#pragma once
#include<iostream>
#include<cmath>

class RefFrame
{
public:
	enum RefFrameType {				//�ο����
		WGS84 = 0,					//GPS
		ITRF,                       //ITRS
		CGCS2000,                   //BDS
		GTRF,                       //Galileo
		PZ90                        //GLONASS
	};
	//-------------field-------------
	RefFrameType    ref;
	double          a;                      //������ unit��m
	double          b;                      //�̰���
	double          f;                      //����
	double          esq;					//��һƫ����ƽ��
	double          GM;                     //������������ unit:m^3/s^2
	double          g_a;                    //����������� unit:m/s^2
	double          U_0;                    //������������λ unit:m^2/s^2
	//-------------constructor-------------
	RefFrame(RefFrameType r = WGS84);
	RefFrame(const RefFrame& right);
	void setFrame(RefFrameType r);
	//-------------function-------------
	double computeR_N(double B);                  //����î��Ȧ�뾶N
	double computeR_M(double B);                  //��������Ȧ�뾶M
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