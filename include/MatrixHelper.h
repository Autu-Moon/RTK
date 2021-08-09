#pragma once
#include<math.h>
#include<iostream>
using namespace std;
class MatrixHelper
{
public://��ģ����ĺ���ֻ��д��ͷ�ļ�
	template<typename T>
	int getLength(T& arr)
	{
		//������鳤��
		return sizeof(arr) / sizeof(arr[0]);
	}

	template<typename T, typename U>
	bool Vector_Plus(const T& vec1, const U& vec2, double Result[])
	{
		/*�����ӷ�*/
		if (getLength(vec1) == getLength(vec2)) {
			//��ӦԪ�����
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
		/*��������*/
		if (getLength(vec1) == getLength(vec2)) {
			//��ӦԪ�����
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
		/*�������*/
		Result = 0;
		if (getLength(vec1) == getLength(vec2)) {
			//��ӦԪ����˺��ۼ�
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
		/*�������*/
		Result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
		Result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
		Result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
		return true;
	}

	template<typename T, typename U>
	bool Matrix_Plus(const T& mat1, int col1, int row1,
		const U& mat2, int col2, int row2, double Result[])
	{
		/*����ӷ�*/

		if (col1 == col2 && row1 == row2) {
			//��ӦԪ�����(����Ҫ�þ�����˼�� ���鼴��)
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
		/*�������*/

		if (col1 == col2 && row1 == row2) {
			//��ӦԪ�����(����Ҫ�þ�����˼�� ���鼴��)
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
		/*����˷�*/
		//�˷�����(���վ�����˼��)
		if (col1 == row2) {
			for (int k = 0; k < row1; k++)
			{
				for (int j = 0; j < col2; j++)
				{
					for (int i = 0; i < col1; i++)
					{
						//i��mat1��Ԫ�ص�λ��,Ҳ��mat2������
						//j��mat2��Ԫ�ص�λ��
						//k��mat1������
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
		/*����ת��*/
		
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
		/*��������*/
		//nΪ����������
		int i, j, k, l, u, v, is[10], js[10];   /* matrix dimension <= 10 */
		double d, p;
		if (n <= 0 || getLength(a) != n * n)
		{
			printf("Error dimension in MatrixInv!\n");
			//exit(EXIT_FAILURE);
			return false;
		}

		/* ���������ֵ���������Result�������Result�������棬a���󲻱� */
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
			for (i = k; i < n; i++)   /* �������½Ƿ�������Ԫ�ص�λ�� */
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

			if (d < DBL_EPSILON)   /* ��Ԫ�ؽӽ���0�����󲻿��� */
			{
				printf("Divided by 0 in MatrixInv!\n");
				return 0;
			}

			if (is[k] != k)  /* ����Ԫ�����ڵ��������½Ƿ�������н��е��� */
			{
				for (j = 0; j < n; j++)
				{
					u = k * n + j;//���½Ƿ�������
					v = is[k] * n + j;//��Ԫ����
					p = Result[u];
					Result[u] = Result[v];
					Result[v] = p;
				}
			}

			if (js[k] != k)  /* ����Ԫ�����ڵ��������½Ƿ�������н��е��� */
			{
				for (i = 0; i < n; i++)
				{
					u = i * n + k;//���½Ƿ�������
					v = i * n + js[k];//��Ԫ����
					p = Result[u];
					Result[u] = Result[v];
					Result[v] = p;
				}
			}

			l = k * n + k;
			Result[l] = 1.0 / Result[l];  /* �����б任 */
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

		for (k = n - 1; k >= 0; k--)  /* ����������е������»ָ� */
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