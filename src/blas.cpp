#include <dymp/blas.h>

#include <algorithm>

namespace dymp{

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////

Vector::Vector(){
	n      = 0;
	vh     = 0;
}

Vector::Vector(int _n){
	n      = 0;
	vh     = 0;
	Allocate(_n);
}

Vector::~Vector(){
	Delete();
}

Vector Vector::SubVector(int ofst, int sz){
	Vector y;
	y.vh = vh + ofst;
	y.n  = sz;
	return y;
}

void Vector::Delete(){
	//if(vh) delete[] vh;
}

void Vector::Allocate(int _n){
	Delete();
	vh = new double[_n];
		
	n = _n;
}

void Vector::Resize(int _n){
	n = _n;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

Matrix::Matrix(){
	m  = 0;
	n  = 0;
	l  = 0;
	vh = 0;
}

Matrix::Matrix(int _m, int _n){
	m  = 0;
	n  = 0;
	l  = 0;
	vh = 0;
	Allocate(_m, _n);
}

Matrix::~Matrix(){
	//if(vh) delete[] vh;
}

void Matrix::Delete(){
	if(vh) delete[] vh;
}

Matrix Matrix::SubMatrix(int row, int col, int _m, int _n){
	Matrix y;
	y.vh = vh + row + l*col;
	y.m = _m;
	y.n = _n;
	y.l =  l;
	return y;
}

Vector Matrix::Col(int col){
	assert(0 <= col && col < n);
	Vector c;
	c.vh = vh + l*col;
	c.n = m;
	return c;
}

void Matrix::Allocate(int _m, int _n){
	Delete();
	vh = new double[_m*_n];
	
	m = _m;
	n = _n;
	l = _m;
}

void Matrix::Resize(int _m, int _n){
	m = _m;
	n = _n;
}

ostream& operator<<(ostream& os, Vector& v){
	for(int i = 0; i < v.n; i++)
		os << v(i) << " ";
	return os;
}

ostream& operator<<(ostream& os, Matrix& m){
	for(int i = 0; i < m.m; i++)for(int j = 0; j < m.n; j++)
		os << m(i,j) << " ";
	return os;
}

}
