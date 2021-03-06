// test
#include "matrix.h"
#include <iostream>

using namespace std;

int main() {
	// test case
	cout << "--------------------------------------------------------------------------------" << endl;
	cout << "TEST----------------------------------------------------------------------------" << endl;
	cout << "--------------------------------------------------------------------------------" << endl;


	int rows = 3;
	int cols = 4;
	double data[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36 };
	double m[] = { 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 };
	double positive[] = { 2, 4, -3, 4, 14, -9, -3, -9, 12 };
	double lower[] = { 1,0,0,0,2,3,0,0,4,5,6,0,7,8,9,10 };
	if (sizeof(data) / sizeof(data[0]) != 6 * 6) {
		cout << "Matrix dimensions do not match!" << endl;
		return 0;
	}
	Matrix a = Matrix(rows, cols, data);
	Matrix b = Matrix(rows, cols, m);
	Matrix c = Matrix(cols, rows, data);
	Matrix l = Matrix(5, 6, data);
	Matrix s = Matrix(3, 3, positive);
	Matrix low = Matrix(4, 4, lower);


	cout << "--------------------a" << endl;
	a.print();
	cout << "--------------------b" << endl;
	b.print();
	cout << "--------------------c" << endl;
	c.print();
	cout << "--------------------s" << endl;
	s.print();
	cout << "--------------------low" << endl;
	low.print();
	cout << "--------------------transpose" << endl;
	Matrix t = a.transpose();
	t.print();
	cout << "--------------------+" << endl;
	Matrix d = a + b;
	d.print();
	cout << "--------------------+-+" << endl;
	Matrix e = a + b + d;
	e.print();
	cout << "--------------------'-'" << endl;
	Matrix f = a - b;
	f.print();
	cout << "--------------------'-'-'-'" << endl;
	Matrix g = a - b - f;
	g.print();
	cout << "--------------------*" << endl;
	Matrix h = a * c;
	h.print();
	cout << "--------------------*-*" << endl;
	Matrix i = a * c * h;
	i.print();
	cout << "--------------------Cholesky Decomposition" << endl;
	Matrix L = s.cholesky_decomposition();
	L.print();
	cout << "--------------------Determinant" << endl;
	double det = s.determinant();
	cout << det << endl;
	cout << "--------------------Identity" << endl;
	Matrix I = Matrix::identity(10);
	I.print();
	cout << "--------------------lower matrix inverse" << endl;
	Matrix Li = low.lower_triangular_inverse();
	Li.print();
	cout << "--------------------inverse" << endl;
	Matrix inv = s.inverse();
	inv.print();


 	cout << "--------------------------------------------------------------------------------" << endl;
	cout << "END-----------------------------------------------------------------------------" << endl;
	cout << "--------------------------------------------------------------------------------" << endl;

}