#include "LeastSquare.h"
bool LeastSquare::Solve(const Matrixd& B, const Matrixd& P, const Matrixd& L,Matrixd& x) {
	try {
		x = (B.Transpose() * P * B).Inverse() * B.Transpose() * P * L;
	}
	catch (Exception e) {
		return false;
	}
	return true;
}
