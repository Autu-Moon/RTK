#pragma once
#include<iostream>
#include<math.h>
#include"MathUtils.h"

struct Hopefield_Element
{
	//���շƶ���ģ�ͱ�׼����Ԫ��
	double H0;						//��ƽ�� m
	double T0;						//�¶�   ����ѧ�¶�
	double P0;						//��ѹ   mbar
	double RH0;						//���ʪ��	
	//Ĭ�ϲ���
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
