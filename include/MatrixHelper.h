#pragma once
#include<math.h>
#include<iostream>
using namespace std;
class MatrixHelper
{
public://用模板类的函数只能写在头文件
	template<typename T>
	int getLength(T& arr)
	{
		//获得数组长度
		return sizeof(arr) / sizeof(arr[0]);
	}

	template<typename T, typename U>
	bool Vector_Plus(const T& vec1, const U& vec2, double Result[])
	{
		/*向量加法*/
		if (getLength(vec1) == getLength(vec2)) {
			//对应元素相加
			for (int i = 0; i < getLength(vec1); i++)
			{
				*(Result + i) = vec1[i] + vec2[i];
			}
			return true;
		}
		else {
			return false;
		}
	}

	template<typename T, typename U>
	bool Vector_Subtract(const T& vec1, const U& vec2, double Result[])
	{
		/*向量减法*/
		if (getLength(vec1) == getLength(vec2)) {
			//对应元素相减
			for (int i = 0; i < getLength(vec1); i++)
			{
				*(Result + i) = vec1[i] - vec2[i];
			}
			return true;
		}
		else {
			return false;
		}
	}

	template<typename T, typename U>
	bool Vector_PointMult(T& vec1, U& vec2, double &Result)
	{
		/*向量点乘*/
		Result = 0;
		if (getLength(vec1) == getLength(vec2)) {
			//对应元素相乘后累加
			for (int i = 0; i < getLength(vec1); i++)
			{
				Result = Result + vec1[i] * vec2[i];
			}
			return true;
		}
		else {
			return false;
		}
	}

	template<typename T, typename U>
	bool Vector_CrossMult(const T vec1[3], const U vec2[3], double Result[])
	{
		/*向量叉乘*/
		Result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
		Result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
		Result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
		return true;
	}

	template<typename T, typename U>
	bool Matrix_Plus(const T& mat1, int col1, int row1,
		const U& mat2, int col2, int row2, double Result[])
	{
		/*矩阵加法*/

		if (col1 == col2 && row1 == row2) {
			//对应元素相减(不需要用矩阵来思考 数组即可)
			for (int i = 0; i < col1*row1; i++)
			{
				Result[i] = mat1[i] + mat2[i];
			}
			return true;
		}
		return false;
	}

	template<typename T, typename U>
	bool Matrix_Subtract(const T& mat1, int col1, int row1,
		const U& mat2, int col2, int row2, double Result[])
	{
		/*矩阵减法*/

		if (col1 == col2 && row1 == row2) {
			//对应元素相减(不需要用矩阵来思考 数组即可)
			for (int i = 0; i < col1*row1; i++)
			{
				Result[i] = mat1[i] - mat2[i];
			}
			return true;
		}
		return false;
	}

	template<typename T, typename U>
	bool Matrix_Mult(const T& mat1, const int &col1,const int &row1,
		const U& mat2, const int &col2,const int &row2, double *Result)
	{
		/*矩阵乘法*/
		//乘法计算(按照矩阵来思考)
		if (col1 == row2) {
			for (int k = 0; k < row1; k++)
			{
				for (int j = 0; j < col2; j++)
				{
					for (int i = 0; i < col1; i++)
					{
						//i是mat1行元素的位次,也是mat2的行数
						//j是mat2列元素的位次
						//k是mat1的行数
						Result[j + k * col2] += mat1[i + k * col1] * mat2[j + i * col2];
					
					}
				}
			}
			return true;
		}
		return false;
	}

	template<typename T>
	bool Matrix_Transpose(const T& mat, int col, int row,double Result[])
	{
		/*矩阵转置*/
		
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
			{
				//Result(i,j)=mat(j,i)
				Result[i + j * row] = mat[j + i * col];
			}
		}
		return true;
	}

	template<typename T>
	bool Matrix_Inv(const T& a,const int n, double Result[])
	{
		/*矩阵求逆*/
		//n为列数或行数
		int i, j, k, l, u, v, is[10], js[10];   /* matrix dimension <= 10 */
		double d, p;
		if (n <= 0 || getLength(a) != n * n)
		{
			printf("Error dimension in MatrixInv!\n");
			//exit(EXIT_FAILURE);
			return false;
		}

		/* 将输入矩阵赋值给输出矩阵Result，下面对Result矩阵求逆，a矩阵不变 */
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				Result[i*n + j] = a[i*n + j];
			}
		}

		for (k = 0; k < n; k++)
		{
			d = 0.0;
			for (i = k; i < n; i++)   /* 查找右下角方阵中主元素的位置 */
			{
				for (j = k; j < n; j++)
				{
					l = n * i + j;
					p = fabs(Result[l]);
					if (p > d)
					{
						d = p;
						is[k] = i;
						js[k] = j;
					}
				}
			}

			if (d < DBL_EPSILON)   /* 主元素接近于0，矩阵不可逆 */
			{
				printf("Divided by 0 in MatrixInv!\n");
				return 0;
			}

			if (is[k] != k)  /* 对主元素所在的行与右下角方阵的首行进行调换 */
			{
				for (j = 0; j < n; j++)
				{
					u = k * n + j;//右下角方阵首行
					v = is[k] * n + j;//主元素行
					p = Result[u];
					Result[u] = Result[v];
					Result[v] = p;
				}
			}

			if (js[k] != k)  /* 对主元素所在的列与右下角方阵的首列进行调换 */
			{
				for (i = 0; i < n; i++)
				{
					u = i * n + k;//右下角方阵首列
					v = i * n + js[k];//主元素列
					p = Result[u];
					Result[u] = Result[v];
					Result[v] = p;
				}
			}

			l = k * n + k;
			Result[l] = 1.0 / Result[l];  /* 初等行变换 */
			for (j = 0; j < n; j++)
			{
				if (j != k)
				{
					u = k * n + j;
					Result[u] = Result[u] * Result[l];
				}
			}
			for (i = 0; i < n; i++)
			{
				if (i != k)
				{
					for (j = 0; j < n; j++)
					{
						if (j != k)
						{
							u = i * n + j;
							Result[u] = Result[u] - Result[i*n + k] * Result[k*n + j];
						}
					}
				}
			}
			for (i = 0; i < n; i++)
			{
				if (i != k)
				{
					u = i * n + k;
					Result[u] = -Result[u] * Result[l];
				}
			}
		}

		for (k = n - 1; k >= 0; k--)  /* 将上面的行列调换重新恢复 */
		{
			if (js[k] != k)
			{
				for (j = 0; j < n; j++)
				{
					u = k * n + j;
					v = js[k] * n + j;
					p = Result[u];
					Result[u] = Result[v];
					Result[v] = p;
				}
			}
			if (is[k] != k)
			{
				for (i = 0; i < n; i++)
				{
					u = i * n + k;
					v = is[k] + i * n;
					p = Result[u];
					Result[u] = Result[v];
					Result[v] = p;
				}
			}
		}
		return true;
	}
};