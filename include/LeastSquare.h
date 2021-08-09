#pragma once
#include"MathUtils.h"
#include"MyMatrix.h"
using namespace MyMatrix;
class LeastSquare
{
public:
	//------------field-------------
	//---------set function---------
	//--------solve function--------
	bool Solve(const Matrixd& B, const Matrixd& P, const Matrixd& L, Matrixd& x);
	//-----------function-----------
	//---------get function---------
	//----------costructor----------
	LeastSquare() = default;
	~LeastSquare() = default;
};

