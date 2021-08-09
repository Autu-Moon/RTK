#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <string>
#include <vector>

namespace
{
    // 辅助代码*
    // =========================================
    // 用来判断类别是否相同
    // =========================================
    // 模板特化
    template <typename T, typename U>
    struct SameType
    {
        static const bool isSame = false;
    };

    template <typename T>
    struct SameType<T, T>
    {
        static const bool isSame = true;
    };

    // 用来判断类别是否是复数
    // =========================================
    template <typename T>
    struct ComplexType {
        static const bool isComplex = false;
    };

    template <typename T>
    struct ComplexType <std::complex<T> > {
        static const bool isComplex = true;
    };

    class Exception
    {
    public:
        explicit Exception(const std::string &_m) : message(_m)
        {
        }

        void printMessage() const
        {
            std::cout << message << std::endl;
        }

    private:
        std::string message;
    };
    // 辅助代码结束
    // =========================================
} 

namespace MyMatrix
{
    template <typename T>
    class Matrix;

    template <typename S>
    std::ostream &operator<<(std::ostream &out, const Matrix<S> &mat);

    template <typename K,typename M>
    Matrix<M> operator*(const K &left, const Matrix<M> &right);
    template <typename K, typename M>
    Matrix<K> operator*(const Matrix<K> &left, const M &right);
    template <typename K, typename M>
    Matrix<M> operator*(const Matrix<K> &left, const Matrix<M> &right);

    template <typename S>
    Matrix<S> operator^(const Matrix<S>& mat, unsigned exponent);


    template <typename T>
    class Matrix
    {
    public:
        //用于构造
        Matrix() = default;
        Matrix(Matrix<T> &&other);
        Matrix(const Matrix<T> &other);
        Matrix(size_t _n);
        Matrix(size_t _x, size_t _y);
        Matrix(const std::vector<std::vector<T> > &dvec);
        Matrix(std::vector<std::vector<T> > &&dvec);


        //工厂函数
        static Matrix<T> Eye(size_t _x, size_t _y);
        static Matrix<T> Eye(size_t _n);
        static Matrix<T> Ones(size_t _x, size_t _y);
        static Matrix<T> Ones(size_t _n);
		static Matrix<T> Zeros(size_t _x, size_t _y);
		static Matrix<T> Zeros(size_t _n);
        //操作矩阵
        static Matrix<T> HStack(const Matrix<T> &left, const Matrix<T> &right);
        static Matrix<T> VStack(const Matrix<T> &top, const Matrix<T> &bottom);
        Matrix<T> Cut(size_t rs, size_t cs, size_t rn, size_t cn) const;
		void SetBlock(size_t rs, size_t cs, size_t rn, size_t cn, const Matrix<T> &block);
        void inline Clear();
        void inline Swap(Matrix<T> &rhs);
        void inline SetLineBreakNum(size_t num);
        void ImportFile(std::string filename);
		void ImportFile(const char *filename);
		void ImportArray(T *array, size_t r, size_t c);
        void ImportArrayMat(T *array, size_t r, size_t c);
        void SetConstant(T var);
        void SetZeros();

        void SetRand(T lower, T upper);
        void SetRand(T lower, T upper, std::mt19937 *_genPtr);
        void SetNormalRand(T mean, T sigma);
        void SetNormalRand(T mean, T sigma, std::mt19937 *_genPtr);

        //用于获取数据
        std::vector<std::vector<T>> &GetData();
        const std::vector<std::vector<T> > &GetData() const;
        T inline GetData(size_t n) const;

        size_t inline Col() const;
        size_t inline Row() const;
        Matrix<double> ToDouble() const;
        void ExportArray(size_t size, T *array) const;

        //用于判断
        bool inline IsEmpty() const;
        bool inline IsSquare() const;
        bool IsSingular() const;

        //向量运算
        static T inline Average(const std::vector<T> &vec);
        static T VecDotProduct(const std::vector<T> lhs, std::vector<T> rhs);

        //矩阵运算
        T inline GetMax() const;
        T inline GetMin() const;
        T inline Average() const;

        Matrix<T> Transpose() const;
        Matrix<double> Inverse() const;
        std::complex<double> Det() const;

        //运算符重载
        std::vector<T> &operator[](size_t index);
        const std::vector<T> &operator[](size_t index) const;
        Matrix<T> operator=(const Matrix<T> &other);
        Matrix<T> operator=(Matrix<T> &&other);
        Matrix<T> operator+(const Matrix<T> &other) const;
        Matrix<T> operator-(const Matrix<T> &other) const;
        Matrix<T> &operator+=(const Matrix<T> &other);
        Matrix<T> &operator-=(const Matrix<T> &other);
        Matrix<T> &operator*=(const Matrix<T> &other);
        Matrix<T> &operator*=(const T &other);
        
        //友元模板实现
        //输出流
        template <typename S>
        friend std::ostream &operator<<(std::ostream &out, const Matrix<S> &mat);

        //operator*
        template <typename K,typename M>
        friend Matrix<M> operator* (const K &left, const Matrix<M> &right);
        template <typename K, typename M>
        friend Matrix<K> operator*(const Matrix<K> &left, const M &right);
        template <typename K,typename M>
        friend Matrix<M> operator* (const Matrix<K> &left, const Matrix<M> &right);

        template <typename S>
        friend Matrix<S> inline operator ^ (const Matrix<S>& mat, unsigned exponent);

    private:
        std::vector<std::vector<T> > data;

        //静态数据声明
        static size_t LineBreakNum;  //输出时的换行数

        //隐藏部分函数
        void Mul(Matrix<T> &ret, const Matrix<T> &other) const;
        void SMul(Matrix<T> &ret, const Matrix<T> &other) const;
        void StrassenMul(size_t rs, size_t re, size_t cs, size_t ce, const Matrix<T> &other, Matrix<T> &ret) const;

        Matrix<double> ConvertToHessenberg() const;
        Matrix<std::complex<double> > GetEig(double eps = 1e-12, unsigned LOOP = 100000) const;
    };

    template <typename T>
    size_t Matrix<T>::LineBreakNum = 1;

    /**
     * 构造函数：移动构造
     */
    template <typename T>
    Matrix<T>::Matrix(Matrix<T> &&other)
    {
        data.swap(other.data);
    }

    /**
     * 构造函数：拷贝构造
     */
    template <typename T>
    Matrix<T>::Matrix(const Matrix<T> &other)
    {
        data = other.GetData();
    }

    /**
     * 构造函数：建立一个n行n列的空矩阵
     */
    template <typename T>
    Matrix<T>::Matrix(size_t _n)
    {
        assert(_n > 0);
        std::vector<std::vector<T> > temp(_n, std::vector<T>(_n));
        data = temp;
    }

    /**
     * 构造函数：根据行列数构造矩阵
     */
    template <typename T>
    Matrix<T>::Matrix(size_t _x, size_t _y)
    {
        assert(_x > 0);
        assert(_y > 0);
        std::vector<std::vector<T> > temp(_x, std::vector<T>(_y));
        data = temp;
    }

    /**
     * 构造函数：调用vector拷贝方法，深拷贝
     */
    template <typename T>
    Matrix<T>::Matrix(const std::vector<std::vector<T> > &dvec)
    {
        assert(dvec.size() > 0);
        assert(dvec[0].size() > 0);
        data = dvec;
    }

    /**
     * 构造函数：调用vector拷贝方法，深拷贝
     */
    template <typename T>
    Matrix<T>::Matrix(std::vector<std::vector<T> > &&dvec)
    {
        assert(dvec.size() > 0);
        assert(dvec[0].size() > 0);
        data = dvec;
    }
    
    

    //从文件输入
    template <typename T>
    void Matrix<T>::ImportFile(std::string filename) 
    {
        assert(std::freopen(filename.c_str(), "r", stdin) != nullptr);

        for (size_t i = 0; i < Row(); ++i) 
        {
            for (size_t j = 0; j < Col(); ++j) 
            {
                std::cin >> data[i][j]; 
            }
        }

        assert(std::freopen("CON", "r", stdin) != nullptr);
    }

    //从文件输入
    template <typename T>
    void Matrix<T>::ImportFile(const char *filename) 
    {
        assert(std::freopen(filename, "r", stdin) != nullptr);

        for (size_t i = 0; i < Row(); ++i) 
        {
            for (size_t j = 0; j < Col(); ++j) 
            {
                std::cin >> data[i][j]; 
            }
        }

        assert(freopen("CON", "r", stdin) != nullptr);
    }
	/*
	 * 导入一个一维数组
	 */
	template<typename T>
	void Matrix<T>::ImportArray(T *array, size_t r, size_t c) {
		for (size_t j = 0; j < c; j++)
		{
			for (size_t i = 0; i < r; i++)
			{
				data[i][j] = array[i + j];
			}
		}
	}

    /**
     * 导入一个一维数组形式的矩阵
	 * 按照r行 c列 的格式
     */
    template <typename T>
    void Matrix<T>::ImportArrayMat(T *array, size_t r, size_t c)
    {
        //assert(!(r > Row()));
        //assert(!(c > Col()));
        for (size_t i = 0; i < r; i++)
        {
            for (size_t j = 0; j < c; j++)
            {
                data[i][j] = array[i * r + j];
            }
        }
    }

    /**
     * 工厂函数：生成一个指定大小的单位阵
     * @param  _x 行数
     * @param  _y  列数
     * @return  一个单位矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::Eye(size_t _x, size_t _y)
    {
        Matrix<T> mat(_x, _y);

        T temp = _x < _y ? _x : _y;
        for (size_t i = 0; i < temp; i++)
        {
            mat.data[i][i] = 1;
        }
        return mat;
    }

    /**
     * 工厂函数：生成一个指定大小的单位阵
     * @param  _n 行数
     * @return  一个单位矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::Eye(size_t _n)
    {
        Matrix<T> mat(_n, _n);

        for (size_t i = 0; i < _n; i++)
        {
            mat.data[i][i] = 1;
        }

        return mat;
    }

    /**
     * 工厂函数：生成一个指定大小的全为1的矩阵
     * @param  _x 行数
     * @param  _y  列数
     * @return  一个全为1的矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::Ones(size_t _x, size_t _y)
    {
        Matrix<T> mat(_x, _y);

        for (size_t i = 0; i < _x; i++)
        {
            for (size_t j = 0; j < _y; j++)
            {
                mat.data[i][j] = 1;
            }
        }

        return mat;
    }

    /**
     * 工厂函数：生成一个指定大小的全为1的矩阵
     * @param  _n 行行数
     * @return  一个全为1的矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::Ones(size_t _n)
    {
        Matrix<T> mat(_n, _n);

        for (size_t i = 0; i < _n; i++)
        {
            for (size_t j = 0; j < _n; j++)
            {
                mat.data[i][j] = 1;
            }
        }

        return mat;
    }
	/**
	 * 工厂函数：生成一个指定大小的全为0的矩阵
	 * @param  _x 行数
	 * @param  _y  列数
	 * @return  一个全为0的矩阵
	 */
	template <typename T>
	Matrix<T> Matrix<T>::Zeros(size_t _x, size_t _y)
	{
		Matrix<T> mat(_x, _y);

		for (size_t i = 0; i < _x; i++)
		{
			for (size_t j = 0; j < _y; j++)
			{
				mat.data[i][j] = 0;
			}
		}

		return mat;
	}

	/**
	 * 工厂函数：生成一个指定大小的全为0的矩阵
	 * @param  _n 行行数
	 * @return  一个全为0的矩阵
	 */
	template <typename T>
	Matrix<T> Matrix<T>::Zeros(size_t _n)
	{
		Matrix<T> mat(_n, _n);

		for (size_t i = 0; i < _n; i++)
		{
			for (size_t j = 0; j < _n; j++)
			{
				mat.data[i][j] = 0;
			}
		}

		return mat;
	}
    /**
     * 水平方向上堆叠
     */
    template <typename T>
    Matrix<T> Matrix<T>::HStack(const Matrix<T> &left,const Matrix<T> &right)
    {
        size_t leftRow = left.Row();
        size_t leftCol = left.Col();
        size_t rightRow = right.Row();
        size_t rightCol = right.Col();

        assert(leftRow == rightRow);

        Matrix<T> mat(leftRow, leftCol + rightCol);
        for (size_t i = 0; i < leftRow; i++)
        {
            for (size_t j = 0; j < leftCol; j++)
            {
                mat[i][j] = left[i][j];
            }
        }
        for (size_t i = 0; i < rightRow; i++)
        {
            for (size_t j = 0; j < rightCol; j++)
            {
                mat[i][leftCol+j] = right[i][j];
            }
        }
        return mat;
    }

    /**
     * 竖直方向上堆叠
     */
    template <typename T>
    Matrix<T> Matrix<T>::VStack(const Matrix<T> &top, const Matrix<T> &bottom)
    {
        size_t topRow = top.Row();
        size_t bottomRow = bottom.Row();

        assert(top.Col() == bottom.Col());
        Matrix<T> mat(topRow + bottomRow, top.Col());

        for (size_t i = 0; i < topRow; i++)
        {
            mat[i] = top[i];
        }
        for (size_t i = 0; i < bottomRow; i++)
        {
            mat[topRow + i] = bottom[i];
        }
        return mat;
        
    }

    /**
     * 从矩阵中取一部分
     * 从 rs,cs开始，大小为rn行 cn列
     * 返回新的另一个矩阵，不影响原矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::Cut(size_t rs, size_t cs, size_t rn, size_t cn) const
    {
        assert(!IsEmpty());
		assert(rs + rn <= Row() && rs >= 0 && cs + cn <= Col() && cs >= 0);
		Matrix<T> ret(rn, cn);

        for (size_t i = rs, ri = 0; ri < rn; ++i, ++ri)
        {
            for (size_t j = cs, rj = 0; rj < cn; ++j, ++rj)
            {
                ret[ri][rj] = data[i][j];
            }
        }
        return ret;
    }
	/**
	 * 为矩阵的一部分赋值
	 * 从 rs,cs开始，长度为rn , cn
	 * 返回新的另一个矩阵，不影响原矩阵
	 */
	template <typename T>
	void Matrix<T>::SetBlock(size_t rs, size_t cs, size_t rn, size_t cn, const Matrix<T> &block) {
		assert(block.Row() <= rn && block.Col() <= cn);
		assert(this->Row() - rs >= rn && this->Col() - cs >= cn);

		for (size_t i = rs, ri = 0; ri < rn; ++i, ++ri)
		{
			for (size_t j = cs, rj = 0; rj < cn; ++j, ++rj)
			{
				this->data[i][j] = block[ri][rj];
			}
		}
	}

    /**
     * 矩阵清空
     * 只删除数据 不清空内存
     */
    template <typename T>
    void Matrix<T>::Clear()
    {
        data.clear();
    }

    /**
     * 调用vector的swap方法，和右端矩阵元素整体调换
     */
    template <typename T>
    void Matrix<T>::Swap(Matrix<T> &rhs)
    {
        data.swap(rhs.GetData());
    }

    /**
    * 控制矩阵输出时的换行数，默认为1
    */
    template <typename T>
    void Matrix<T>::SetLineBreakNum(size_t num)
    {
        Matrix<T>::LineBreakNum = num;
    }

    /**
     * 矩阵清零
     */
    template <typename T>
    void Matrix<T>::SetZeros()
    {
        size_t n = Row();
        size_t m = Col();

        for (size_t i = 0; i < n; i++)
        {
            std::memset(&data[i][0], 0, sizeof(T) * m);
        }
    }

    /**
    * 矩阵设为常数
    */
    template <typename T>
    void Matrix<T>::SetConstant(T var)
    {
        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                data[i][j] = var;
            }
        }
    }

   

    /**
    * 矩阵设为随机数
    */
    template <typename T>
    void Matrix<T>::SetRand(T lower, T upper)
    {
        std::uniform_real_distribution<T> unifDis(lower, upper);
        std::mt19937 *_genPtr = new std::mt19937(std::time(NULL));
        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                data[i][j] = unifDis(*_genPtr);
            }
        }
        delete _genPtr;
    }

    /**
    * 矩阵设为随机数
    */
    template <typename T>
    void Matrix<T>::SetRand(T lower, T upper, std::mt19937 *_genPtr)
    {
        assert(_genPtr != nullptr);

        std::uniform_real_distribution<T> unifDis(lower, upper);
        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                data[i][j] = unifDis(*_genPtr);
            }
        }
    }

    /**
    * 矩阵设为正态分布随机数
    */
    template <typename T>
    void Matrix<T>::SetNormalRand(T mean, T sigma)
    {
        std::mt19937 *_genPtr = new std::mt19937(std::time(NULL));
        std::normal_distribution<T> normDis(mean, sigma);
        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                data[i][j] = normDis(*_genPtr);
            }
        }
        delete _genPtr;
    }

    /**
    * 矩阵设为正态分布随机数
    */
    template <typename T>
    void Matrix<T>::SetNormalRand(T mean, T sigma, std::mt19937 *_genPtr)
    {
        assert(_genPtr != nullptr);

        std::normal_distribution<T> normDis(mean, sigma);
        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                data[i][j] = normDis(*_genPtr);
            }
        }
    }

    /**
     * 获取 data 成员，可用于整块更新
     */
    template <typename T>
    std::vector<std::vector<T> > &Matrix<T>::GetData()
    {
        return data;
    }
    

    /**
     * 以 const 方式获取成员，可用于安全读
     */
    template <typename T>
    const std::vector<std::vector<T> > &Matrix<T>::GetData() const
    {
        return data;
    }

    /**
    * 返回第n个数
    */
    template <typename T>
    T Matrix<T>::GetData(size_t n) const
    {
#ifndef NUMCHECK
#define NUMCHECK
        assert(data.size());
#endif
        assert((Col() * Row()) >= n);
        if (n % Col() == 0)
        {
            return data[n / Col() - 1][Col() - 1];
        }
        return data[n / Col()][n % Col() - 1];
    }

    /**
     * 矩阵的行数
     */
    template <typename T>
    size_t Matrix<T>::Row() const
    {
        return data.size();
    }

    /**
     * 矩阵的列数
     */
    template <typename T>
    size_t Matrix<T>::Col() const
    {
        if (data.size())
        {
            return data[0].size();
        }
        else
        {
            return 0;
        }
    }

    /**
     * 转换为 double 矩阵
     */
    template <typename T>
    Matrix<double> Matrix<T>::ToDouble() const
    {
        Matrix<double> mat(Row(), Col());

        // 未实现：如果矩阵是复数，则只取实数部分，忽略虚数部分
        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                mat[i][j] = static_cast<double>(data[i][j]);
            }
        }

        return mat;
    }
   
    /**
    * 提取矩阵前n个元素到数组
    */
    template <typename T>
    void Matrix<T>::ExportArray(size_t size, T *array) const
    {
        assert(!IsEmpty());
        assert((Col() * Row()) >= size);

        size_t k = 0;
        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                array[k] = data[i][j];
                k++;
                if (k == size)
                {
                    break;
                }
            }
            if (k == size)
            {
                break;
            }
        }
    }

    /**
     * 判断是否是空矩阵
     * @return 1 : 0 -> 空矩阵 : 不是空矩阵
     */
    template <typename T>
    bool Matrix<T>::IsEmpty() const
    {
        return !data.size();
    }

    /**
     * 判断矩阵是否是方阵
     * @return 1 : 0 -> 方阵 : 不是方阵
     */
    template <typename T>
    bool Matrix<T>::IsSquare() const
    {
        return IsEmpty() ? 0 : data.size() == data[0].size();
    }

    /**
     * 判断当前矩阵是否是奇异矩阵
     * @return 1: 是奇异矩阵 0: 不是奇异矩阵
     */
    template <typename T>
    bool Matrix<T>::IsSingular() const 
    {
        std::complex<double> temp = Det();
        
        return std::fabs(temp.real()) < 1e-10 ? true : false;
    }


    /**
     * 一个向量的均值
     */
    template <typename T>
    T Matrix<T>::Average(const std::vector<T> &vec)
    {
        T sum = 0;

        for (T var : vec)
        {
            sum += var;
        }

        return sum / vec.size();
    }

    /**
     * 静态函数：获取两个向量的点乘结果
     * @param  lhs 向量1
     * @param  rhs 向量2
     * @return     double:点乘的结果
     */
    template <typename T>
    T Matrix<T>::VecDotProduct(const std::vector<T> lhs, std::vector<T> rhs)
    {
        T ans = 0;

        for (decltype(lhs.size()) i = 0; i != lhs.size(); i++)
        {
            ans += lhs[i] * rhs[i];
        }

        return ans;
    }

    /**
     * 返回矩阵中最大的元素
     */
    template <typename T>
    T Matrix<T>::GetMax() const
    {
        if (!data.size())
        {
            return static_cast<T>(0);
        }

        T maxv = data[0][0];

        for (size_t i = 0; i < data.size(); i++)
        {
            for (size_t j = 0; j < data[0].size(); j++)
            {
                maxv = data[i][j] < maxv ? maxv : data[i][j];
            }
        }

        return maxv;
    }

    /**
     * 返回矩阵中最小的元素
     */
    template <typename T>
    T Matrix<T>::GetMin() const
    {
        if (IsEmpty())
        {
            return static_cast<T>(0);
        }

        T minv = data[0][0];

        for (size_t i = 0; i < data.size(); i++)
        {
            for (size_t j = 0; j < data[0].size(); j++)
            {
                minv = data[i][j] < minv ? data[i][j] : minv;
            }
        }

        return minv;
    }

    /**
     * 矩阵的均值
     */
    template <typename T>
    T Matrix<T>::Average() const
    {
        if (IsEmpty())
        {
            return static_cast<T>(0);
        }

        T sum = 0;

        for (size_t i = 0; i < data.size(); i++)
        {
            for (size_t j = 0; j < data[0].size(); j++)
            {
                sum += data[i][j];
            }
        }

        return sum / (Row() * Col());
    }

    /**
     * 获得矩阵的转置
     * @return 新的矩阵，内容为原矩阵的转置
     */
    template <typename T>
    Matrix<T> Matrix<T>::Transpose() const
    {
        decltype(data.size()) sizeRow = data.size();

        if (sizeRow == 0)
        {
            std::cerr << "error** Matrix<T>::Transposition -> Empty Matrix!" << std::endl;
        }

        using size_t = decltype(data.size());
        size_t sizeCol = data[0].size();

        Matrix tran(sizeCol, sizeRow);

        for (size_t i = 0; i < sizeRow; i++)
        {
            for (size_t j = 0; j < sizeCol; j++)
            {
                tran.data[j][i] = data[i][j];
            }
        }

        return tran;
    }

    /**
     * 高斯约当法求逆
     */
    template <typename T>
    Matrix<double> Matrix<T>::Inverse() const
    {
        assert(!IsEmpty());
        assert(IsSquare());
        //assert(!IsSingular());
        //throw (Exception("ERROR** Matrix::Inverse -> there is no inverse matrix!"));

        size_t i, j, k, len = Row();
        double maxVal, temp;
        //将A矩阵存放在临时矩阵中
        Matrix<double> TMat;

        if (SameType<T, double>::isSame)
        {
            TMat = *this;
        }
        else
        {
            TMat = this->ToDouble();
        }

        //初始化ans矩阵为单位阵
        Matrix<double> ans = Matrix<double>::Eye(Row(), Col());

        for (i = 0; i < len; i++)
        {
            //寻找主元
            maxVal = TMat[i][i];
            k = i;

            for (j = i + 1; j < len; j++)
            {
                if (std::abs(TMat[j][i]) > std::abs(maxVal))
                {
                    maxVal = TMat[j][i];
                    k = j;
                }
            }

            //如果主元所在行不是第i行，进行行交换
            if (k != i)
            {
                TMat[i].swap(TMat[k]);
                ans[i].swap(ans[k]);
            }

            //消去A的第i列除去i行以外的各行元素
            temp = TMat[i][i];

            for (j = 0; j < len; j++)
            {
                TMat[i][j] = TMat[i][j] / temp; //主对角线上的元素变为1
                ans[i][j] = ans[i][j] / temp;   //伴随计算
            }

            // 遍历行
            for (j = 0; j < len; j++)
            {
                // 不是第i行
                if (j != i)
                {
                    temp = TMat[j][i];

                    // 第j行元素 - i行元素 * j列i行元素
                    for (k = 0; k < len; k++)
                    {
                        TMat[j][k] -= TMat[i][k] * temp;
                        ans[j][k] -= ans[i][k] * temp;
                    }
                }
            }
        }
        return ans;
    }

    /**
     * 求矩阵行列式
     * @return double: 行列式的值
     */
    template <typename T>
    std::complex<double> Matrix<T>::Det() const 
    {
        // 所有特征根的乘积
        auto e = GetEig();
        std::complex<double> ret = e[0][0];

        for (size_t i = 1; i < Col(); ++i)
        {
            ret *= e[i][0];
        }

        return ret;
    }

    
    /**
     * 取矩阵中的某一行
     */
    template <typename T>
    std::vector<T> &Matrix<T>::operator[](size_t index)
    {
        assert(index >= 0);
        assert(index < data.size());
        return data[index];
    }

    /**
     * 常对象取矩阵中的某一行
     */
    template <typename T>
    const std::vector<T> &Matrix<T>::operator[](size_t index) const
    {
        assert(index >= 0);
        assert(index < data.size());
        return data[index];
    }

    /**
     * 拷贝矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator=(const Matrix<T> &other)
    {
        data = other.GetData();
        return *this;
    }

    /**
     * 移动矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator=(Matrix<T> &&other)
    {
        data.swap(other.GetData());
        return *this;
    }

    /**
     * 矩阵相加
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const
    {
        assert(!IsEmpty());
        assert(Row() == other.Row() && Col() == other.Col());

        Matrix<T> ret(Row(), Col());

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                ret[i][j] = data[i][j] + other[i][j];
            }
        }

        return ret;
    }

    /**
     * 矩阵相加赋值
     */
    template <typename T>
    Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &other)
    {
        assert(!IsEmpty());
        assert(Row() == other.Row() && Col() == other.Col());

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                data[i][j] += other[i][j];
            }
        }

        return *this;
    }

    /**
     * 矩阵相减
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const
    {
        assert(!IsEmpty());
        assert(Row() == other.Row() && Col() == other.Col());

        Matrix<T> ret(Row(), Col());

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                ret[i][j] = data[i][j] - other[i][j];
            }
        }

        return ret;
    }

    /**
     * 矩阵相减赋值
     */
    template <typename T>
    Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &other)
    {
        assert(!IsEmpty());
        assert(Row() == other.Row() && Col() == other.Col());

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t j = 0; j < Col(); j++)
            {
                data[i][j] -= other[i][j];
            }
        }

        return *this;
    }

    /**
     * 矩阵 *=
     */
    template <typename T>
    Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &other)
    {
        assert(!IsEmpty());
        assert(Col() == other.Row());
        //Matrix<T> ret(Row(), other.Col());

        // 对称矩阵使用 斯特拉森乘法
        if ((Row() == Col()) && (other.Row() == other.Col()))
        {
            SMul(*this, other);
        }
        else
        {
            // 普通乘法
            Mul(*this, other);
        }

        return *this;
    }

    /**
     *  矩阵 *=
     *  立即数
     */
    template <typename T>
    Matrix<T> &Matrix<T>:: operator*=(const T &other) 
    {
        size_t matRow = Row();
        size_t matCol = Col();
        //Matrix<T> ret(matRow,matCol);

        for (size_t i = 0; i < matRow; i++)
        {
            for (size_t j = 0; j < matCol; j++)
            {
                data[i][j] = other * data[i][j];
            }
            
        }

        return *this;
    }

    /**
     * 友元模板实现输出
     */
    template <typename S>
    std::ostream &operator<<(std::ostream &out, const Matrix<S> &mat)
    {
        //判断要输出的数是否是 double 类型
        bool dFlag = SameType<S, double>::isSame;

        if (dFlag)
        {
            using std::ios;
            using std::setprecision;
            out << setiosflags(ios::right) << setprecision(5);
        }

        for (size_t i = 0; i != mat.data.size(); i++)
        {
            for (size_t j = 0; j != mat.data[i].size(); j++)
            {
                if (dFlag)
                {
                    out << std::setw(12) << mat.data[i][j] << ' ';
                }
                else
                {
                    out << mat.data[i][j] << ' ';
                }
            }

            if (i < mat.data.size() - 1)
            {
                out << '\n';
            }
        }
        out << std::endl;
        for (size_t i = 0; i < Matrix<S>::LineBreakNum; i++)
        {
            out << std::endl;
        }
        return out;
    }

    template <typename K,typename M>
    Matrix<M> operator* (const K &left, const Matrix<M> &right)
    {
        size_t matRow = right.Row();
        size_t matCol = right.Col();
        Matrix<M> ret(matRow,matCol);

        for (size_t i = 0; i < matRow; i++)
        {
            for (size_t j = 0; j < matCol; j++)
            {
                ret[i][j] = (left * right[i][j]);
            }
            
        }

        return ret;
    }

    template <typename K,typename M>
    Matrix<K> operator* (const Matrix<K> &left, const M &right)
    {
        size_t matRow = left.Row();
        size_t matCol = left.Col();
        Matrix<K> ret(matRow, matCol);

        for (size_t i = 0; i < matRow; i++)
        {
            for (size_t j = 0; j < matCol; j++)
            {
                ret[i][j] = (right * left[i][j]);
            }
            
        }

        return ret;
    }

    template <typename K,typename M>
    Matrix<M> operator* (const Matrix<K> &left, const Matrix<M> &right)
    {
        
        assert(left.Col() == right.Row());

        Matrix<K> ret(left.Row(), right.Col());

        // 对称矩阵使用 斯特拉森乘法
        if ((left.Row() ==left. Col()) && (right.Row() == right.Col()))
        {
            left.SMul(ret, right);
        }
        // 普通乘法
        else
        {
            left.Mul(ret, right);
        }
        return ret;
    }

    /**
    * 矩阵快速幂
    */
    template <typename S>
    Matrix<S> operator ^ (const Matrix<S>& mat, unsigned exponent)
    {
        assert(mat.Row() == mat.Col());
        Matrix<S> ans = Matrix<S>::Eye(mat.Row(), mat.Col());
        Matrix<S> src = mat;

        for (; exponent; exponent >>= 1)
        {
            if (exponent & 1)
            {
                ans *= src;
            }
            src *= src;
        }

        return ans;
    }

    /**
     * 矩阵相乘
     * 普通算法
     */
    template <typename T>
    void Matrix<T>::Mul(Matrix<T> &ret, const Matrix<T> &other) const
    {
        assert(!IsEmpty());
        assert(Col() == other.Row());

        for (size_t i = 0; i < Row(); i++)
        {
            for (size_t k = 0; k < Col(); k++)
            {
                for (size_t j = 0; j < other.Col(); j++)
                {
                    ret[i][j] += (data[i][k] * other[k][j]);
                }
            }
        }

        return;
    }

    /**
     * 斯特拉森乘法主函数，两个 n * n 矩阵
     */
    template <typename T>
    void Matrix<T>::SMul(Matrix<T> &ret, const Matrix<T> &other) const
    {
        assert(!IsEmpty());
        assert(Col() == other.Row());

        assert(IsSquare());
        assert(other.IsSquare());
        size_t n = Row();
        StrassenMul(0, n, 0, n, other, ret);
    }

    
    /**
     * 斯特拉森乘法
     * 时间复杂度 n ^ 2.80，性能更好，但参数复杂
     */
    template <typename T>
    void Matrix<T>::StrassenMul(size_t rs, size_t re, size_t cs, size_t ce, const Matrix<T> &other, Matrix<T> &ret) const
    {
        assert(!IsEmpty());
        assert(Col() == other.Row());

        if (re - rs == 2 && ce - cs == 2)
        {
            size_t rs1 = rs + 1;
            size_t cs1 = cs + 1;
            T P1 = data[rs][cs] * (other[rs][cs1] - other[rs1][cs1]);
            T P2 = (data[rs][cs] + data[rs][cs1]) * other[rs1][cs1];
            T P3 = (data[rs1][cs] + data[rs1][cs1]) * other[rs][cs];
            T P4 = data[rs1][cs1] * (other[rs1][cs] - other[rs][cs]);
            T P5 = (data[rs][cs] + data[rs1][cs1]) * (other[rs][cs] + other[rs1][cs1]);
            T P6 = (data[rs][cs1] - data[rs1][cs1]) * (other[rs1][cs] + other[rs1][cs1]);
            T P7 = (data[rs][cs] - data[rs1][cs]) * (other[rs][cs] + other[rs][cs1]);
            ret[rs][cs] = P5 + P4 - P2 + P6;
            ret[rs][cs1] = P1 + P2;
            ret[rs1][cs] = P3 + P4;
            ret[rs1][cs1] = P1 + P5 - P3 - P7;
        }
        else if (re - rs < 2 || rs - rs < 2)
        {
            for (size_t i = rs; i < re; i++)
            {
                for (size_t k = cs; k < ce; k++)
                {
                    for (size_t j = cs; j < ce; j++)
                    {
                        ret[i][j] += data[i][k] * other[k][j];
                    }
                }
            }
        }
        else
        {
            size_t rm = rs + ((re - rs) / 2);
            size_t cm = cs + ((ce - cs) / 2);
            StrassenMul(rs, rm, cs, cm, other, ret);
            StrassenMul(rm, re, cs, cm, other, ret);
            StrassenMul(rs, rm, cm, ce, other, ret);
            StrassenMul(rm, re, cm, ce, other, ret);
        }
    }

    /**
     * 生成当前矩阵的 Hessenberg 形式，以新矩阵返回
     */
    template <typename T>
    Matrix<double> Matrix<T>::ConvertToHessenberg() const 
    {
        assert(!IsEmpty());
        assert(IsSquare());

        Matrix<double> A = ToDouble();

        size_t n = data.size();
        size_t i, j, k;
        Matrix<double> ret(n, n);
        T temp = 0;
        size_t max;

        for (k = 1; k < n - 1; ++k) {
            i = k - 1;
            max = k;
            temp = std::abs(A[k][i]);

            for (j = k + 1; j < n; ++j) {
                if (temp < std::abs(A[j][i])) {
                    temp = std::abs(A[j][i]);
                    max = j;
                }
            }

            ret[0][0] = A[max][i];
            i = max;

            if (ret[0][0]) {
                if (i != k) {
                    for (j = k - 1; j < n; ++j) {
                        temp = A[i][j];
                        A[i][j] = A[k][j];
                        A[k][j] = temp;
                    }

                    for (j = 0; j < n; ++j) {
                        temp = A[j][i];
                        A[j][i] = A[j][k];
                        A[j][k] = temp;
                    }
                }

                for (i = k + 1; i < n; ++i) {
                    temp = A[i][k - 1] / ret[0][0];
                    A[i][k - 1] = 0;

                    for (size_t j = k; j < n; ++j) {
                        A[i][j] -= temp * A[k][j];
                    }

                    for (j = 0; j < n; ++j) {
                        A[j][k] += temp * A[j][i];
                    }
                }
            }
        }
        return A;
    }

    /**
     * 返回矩阵的全部特征根，以复数表示
     * QR 分解法
     */
    template<typename T>
    Matrix<std::complex<double> > Matrix<T>::GetEig(double eps, unsigned LOOP) const
    {
        assert(!IsEmpty());
        assert(IsSquare());

        unsigned loop = LOOP;
        const size_t n = data.size();
        size_t m = n;
        Matrix<double> A = ConvertToHessenberg();
        Matrix<double> ret(n, 2);
        size_t i, j, k, t;
        double tempVar, indexVar, sign, p, q;
        double r, x, s, e, f, z, y, temp;
        double num;

        while (m != 0) {
            t = m - 1;

            while (t > 0) {
                temp = std::abs(A[t - 1][t - 1]);
                temp += std::abs(A[t][t]);
                temp *= eps;

                if (std::abs(A[t][t - 1]) > temp) {
                    --t;
                }
                else {
                    break;
                }
            }

            if (t == m - 1) {
                ret[m - 1][0] = A[m - 1][m - 1];
                ret[m - 1][1] = 0;
                m -= 1;
                loop = LOOP;
            }
            else if (t == m - 2) {
                tempVar = -A[m - 1][m - 1] - A[m - 2][m - 2];
                indexVar = A[m - 1][m - 1] * A[m - 2][m - 2]
                    - A[m - 1][m - 2] * A[m - 2][m - 1];
                num = tempVar * tempVar - 4 * indexVar;
                y = std::sqrt(std::abs(num));

                if (num > 0) {
                    sign = 1;

                    if (tempVar < 0) {
                        sign = -1;
                    }

                    ret[m - 1][0] = -(tempVar + sign * y) / 2;
                    ret[m - 1][1] = 0;
                    ret[m - 2][0] = indexVar / ret[m - 1][0];
                    ret[m - 2][1] = 0;
                }
                else {
                    ret[m - 1][0] = -tempVar / 2;
                    ret[m - 2][0] = -tempVar / 2;
                    ret[m - 1][1] = y / 2;
                    ret[m - 2][1] = -y / 2;
                }

                m -= 2;
                loop = LOOP;
            }
            else {
                if (loop < 1) {
                    return Matrix<std::complex<double> >();
                }

                --loop;
                j = t + 2;

                while (j < m) {
                    A[j][j - 2] = 0;
                    ++j;
                }

                j = t + 3;

                while (j < m) {
                    A[j][j - 3] = 0;
                    ++j;
                }

                k = t;

                while (k < m - 1) {
                    if (k != t) {
                        p = A[k][k - 1];
                        q = A[k + 1][k - 1];

                        if (k != m - 2) {
                            r = A[k + 2][k - 1];
                        }
                        else {
                            r = 0;
                        }
                    }
                    else {
                        tempVar = A[m - 1][m - 1];
                        indexVar = A[m - 2][m - 2];
                        x = tempVar + indexVar;
                        y = tempVar * indexVar - A[m - 2][m - 1] * A[m - 1][m - 2];
                        p = A[t][t] * (A[t][t] - x) + A[t][t + 1] * A[t + 1][t] + y;
                        q = A[t + 1][t] * (A[t][t] + A[t + 1][t + 1] - x);
                        r = A[t + 1][t] * A[t + 2][t + 1];
                    }

                    if (p != 0 || q != 0 || r != 0) {
                        if (p < 0) {
                            sign = -1;
                        }
                        else {
                            sign = 1;
                        }

                        s = sign * std::sqrt(p * p + q * q + r * r);

                        if (k != t) {
                            A[k][k - 1] = -s;
                        }

                        e = -q / s;
                        f = -r / s;
                        x = -p / s;
                        y = -x - f * r / (p + s);
                        num = e * r / (p + s);
                        z = -x - e * q / (p + s);

                        for (j = k; j < m; ++j) {
                            tempVar = A[k][j];
                            indexVar = A[k + 1][j];
                            p = x * tempVar + e * indexVar;
                            q = e * tempVar + y * indexVar;
                            r = f * tempVar + num * indexVar;

                            if (k != m - 2) {
                                tempVar = A[k + 2][j];
                                p += f * tempVar;
                                q += num * tempVar;
                                r += z * tempVar;
                                A[k + 2][j] = r;
                            }

                            A[k + 1][j] = q;
                            A[k][j] = p;
                        }

                        j = k + 3;

                        if (j > m - 2) {
                            j = m - 1;
                        }

                        for (i = t; i < j + 1; ++i) {
                            tempVar = A[i][k];
                            indexVar = A[i][k + 1];
                            p = x * tempVar + e * indexVar;
                            q = e * tempVar + y * indexVar;
                            r = f * tempVar + num * indexVar;

                            if (k != m - 2) {
                                tempVar = A[i][k + 2];
                                p += f * tempVar;
                                q += num * tempVar;
                                r += z * tempVar;
                                A[i][k + 2] = r;
                            }

                            A[i][k + 1] = q;
                            A[i][k] = p;
                        }
                    }

                    ++k;
                }
            }
        }

        // 返回一个复数
        Matrix<std::complex<double> > res(n, 1);

        for (size_t i = 0; i < ret.Row(); ++i) {
            // 判断是否是复数类型
            bool flag = ComplexType<T>::isComplex;

            if (flag) 
            {
                res[i][0] = std::complex<double>(static_cast<std::complex<T>>(ret[i][0]).real(), static_cast<std::complex<T>>(ret[i][1]).real());
            }
            else 
            {
                res[i][0] = std::complex<double>(ret[i][0], ret[i][1]);
            }
        }

        return res;
    }

    using Matrixd = Matrix<double>;
    using Matrixf = Matrix<float>;
    using Matrixi = Matrix<int>;
} 

namespace mym = MyMatrix;

/******************************************************************************************
 * Matrix in C++
 * Computer Science & Technology, WuHan University
 * CopyRight(©) 2020 sinecera. All rights reserved
******************************************************************************************/
