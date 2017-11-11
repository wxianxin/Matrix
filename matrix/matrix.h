// © 2017 XIANXIN WANG ALL RIGHTS RESERVED
// Steven Wang 2017
// c++ 11

#ifndef MATRIX_H_
#define MATRIX_H_

class Matrix {
public:
	Matrix(int rows, int cols);
	Matrix(int rows, int cols, double data[]);
	static Matrix identity(int n);
	void print();
	Matrix transpose();
	Matrix operator+ (const Matrix& b);
	Matrix operator- (const Matrix& b);
	Matrix operator* (const Matrix& b);
	Matrix cholesky_decomposition();
	double* forward_substitution(Matrix L, double* b);
	Matrix lower_triangular_inverse();
	double determinant();
	Matrix inverse();

private:
	int rows_;
	int cols_;
	double* p_;	// pointer to the start of array memory
	double** ptr_;	// pointer to the array of pointers, which point to the start of each row in the matrix
};

#endif