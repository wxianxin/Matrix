// © 2017 XIANXIN WANG ALL RIGHTS RESERVED
// Steven Wang 2017
// c++ 11

#include "matrix.h"
#include <cmath>
#include <iostream>
#include <iomanip>

// To build dynamic array of arrays 
// 1: use vector library
// 2: array of pointers to arrays
// 3: allocate your whole matrix as a single dynamic array, then use (slightly) clever indexing math of your own to access cells.
// contiguous memory allocation + multiple column pointers -> make development much easier and less prone to mistakes

Matrix::Matrix(int rows, int cols)
{
	// check input
	//if (sizeof(data) / sizeof(data[0]) != rows * cols) {
	//	std::cout << sizeof(data) << "Matrix dimensions do not match!" << sizeof(data[0]) << std::endl;
	//}
	this->rows_ = rows;
	this->cols_ = cols;
	this->p_ = new double[rows_ * cols_];
	this->ptr_ = new double*[rows];
	int j = 0;
	for (int i = 0; i < this->rows_ * this->cols_; i++)
	{
		if (i%cols == 0)
		{
			ptr_[j++] = &p_[i];
		}
		p_[i] = 0.0;
	}
}

Matrix::Matrix(int rows, int cols, double data[])
{
	// check input
	//if (sizeof(data) / sizeof(data[0]) != rows * cols) {
	//	std::cout << sizeof(data) << "Matrix dimensions do not match!" << sizeof(data[0]) << std::endl;
	//}
	this->rows_ = rows;
	this->cols_ = cols;
	this->p_ = new double[rows_ * cols_];
	this->ptr_ = new double*[rows];
	int j = 0;
	for (int i = 0; i < this->rows_ * this->cols_; i++)
	{
		if (i%cols == 0)
		{
			ptr_[j++] = &p_[i];
		}
		p_[i] = data[i];
	}
}

// TODO clean memory
// https://stackoverflow.com/questions/23566114/reason-for-double-free-or-corruption
//~Matrix() {
//	if(this->p_){
//		delete[] p_;
//	}
//}

Matrix Matrix::identity(int n)
{
	Matrix I(n, n);
	for (int i = 0; i < n; i++)
	{
		I.ptr_[i][i] = 1;
	}

	return I;
}

void Matrix::print() {
	for (int i = 0; i < rows_; i++) {
		for (int j = 0; j < cols_; j++) {
			std::cout << std::setw(8) << std::setprecision(2) << ptr_[i][j];
		}
		std::cout << std::endl;
	}
}

Matrix Matrix::transpose()
{
	Matrix t(this->cols_, this->rows_);
	for (int i = 0; i < this->rows_; i++) {
		for (int j = 0; j < this->cols_; j++) {
			t.ptr_[j][i] = this->ptr_[i][j];
		}
	}

	return t;
}

Matrix Matrix::operator+(const Matrix& b)
{
	Matrix result = Matrix(this->rows_, this->cols_);
	int j = 0;
	for (int i = 0; i < this->rows_ * this->cols_; i++)
	{
		result.p_[i] = this->p_[i] + b.p_[i];
	}

	return result;
}

Matrix Matrix::operator-(const Matrix& b)
{
	Matrix result = Matrix(this->rows_, this->cols_);
	int j = 0;
	for (int i = 0; i < this->rows_ * this->cols_; i++)
	{
		result.p_[i] = this->p_[i] - b.p_[i];
	}
	return result;
}

Matrix Matrix::operator*(const Matrix& b)
{
	Matrix result = Matrix(this->rows_, b.cols_);
	for (int i = 0; i < this->rows_; i++)
	{
		for (int j = 0; j < b.cols_; j++)
		{
			result.ptr_[i][j] = 0.0;
			for (int c = 0; c < this->cols_; c++) {
				result.ptr_[i][j] += (this->ptr_[i][c] * b.ptr_[c][j]);
			}
		}
	}

	return result;
}

Matrix Matrix::cholesky_decomposition()
{
	// cholesky_decomposition (LL^T Decomposition)
	// Cholesky–Banachiewicz algorithm, about O(n^3/3)
	if (this->rows_ != this->cols_)
	{
		std::cout << "Matrix is not a square matrix!!!" << std::endl;
		// TODO: end function
	}

	Matrix L(this->rows_, this->rows_);
	double sum_i, sum_j;
	for (int i = 0; i < this->rows_; i++)
	{
		sum_i = 0;
		for (int j = 0; j < i; j++)
		{
			sum_j = 0;
			for (int k = 0; k < j; k++)
			{
				sum_j += L.ptr_[i][k] * L.ptr_[j][k];
			}
			L.ptr_[i][j] = (this->ptr_[i][j] - sum_j) / L.ptr_[j][j];

			sum_i += L.ptr_[i][j] * L.ptr_[i][j];
		}
		L.ptr_[i][i] = sqrt(this->ptr_[i][i] - sum_i);
	}

	return L;
}

double* Matrix::forward_substitution(Matrix L, double* b) {
	// https://mobiusfunction.wordpress.com/2010/12/21/the-inverse-of-a-triangular-matrix-using-forward-substitution/
	// https://en.wikipedia.org/wiki/Triangular_matrix#Forward_substitution
	//if (L.rows_ != sizeof(b[]))
	//{
	//}
	double* x = new double[L.rows_];
	double sum;
	for (int i = 0; i < L.rows_; i++)
	{
		sum = 0.0;
		for (int c = 0; c < i; c++)
		{
			sum += L.ptr_[i][c] * x[c];
		}

		x[i] = (b[i] - sum) / L.ptr_[i][i];
	}

	return x;
}

Matrix Matrix::lower_triangular_inverse() {
	// For convenience, define: L" = L^(-1)
	// 0. N-by-N: L * L" = I
	// 1. L * L"[ ,i] = I[ ,i] : We get a system of N equations because there is n cols in I
	// 2. For each equation, use forward substitution

	Matrix Li(this->cols_, this->cols_);
	Matrix I = identity(this->cols_);

	// this new is useless
	double* x;
	double* b = new double[this->cols_];
	for (int i = 0; i < this->cols_; i++)
	{
		for (int c = 0; c < this->cols_; c++)
		{
			b[c] = I.ptr_[c][i];
		}
		x = forward_substitution(*this, b);
		for (int c = 0; c < this->cols_; c++)
		{
			Li.ptr_[c][i] = x[c];
		}

		delete[] x;
	}

	delete[] b;

	return Li;
}

double Matrix::determinant()
{
	// calculate determinant
	double det = 1.0;
	Matrix L = this->cholesky_decomposition();
	for (int i = 0; i < L.rows_; i++)
	{
		det *= L.ptr_[i][i] * L.ptr_[i][i];
	}

	return det;
}

Matrix Matrix::inverse()
	{
		// https://scicomp.stackexchange.com/questions/19628/what-is-the-fastest-algorithm-for-computing-the-inverse-matrix-and-its-determina
		// explicit matrix inverse
		// 1. Cholesky Decomposition
		// 2. inverse lower triangular matrix

		Matrix L = this->cholesky_decomposition();
		Matrix Li = L.lower_triangular_inverse();
		Matrix inv = Li.transpose() * Li;

		return inv;
	}
