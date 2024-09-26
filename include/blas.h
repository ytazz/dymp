#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <assert.h>

#ifdef _WIN32
# include <mkl.h>
#else
# include <mkl.h>
#endif

#include <Eigen/Eigen>
//using namespace Eigen;

namespace dymp{;

struct Vector{
	int     n;
	double* vh;
	
	void Delete  ();
	void Allocate(int _n);
	void Resize  (int _n);

	double& operator()(int i)     { return vh[i]; }
	double  operator()(int i)const{ return vh[i]; }

	Vector SubVector(int ofst, int sz);

	 Vector();
	 Vector(int _n);
	~Vector();
};

struct Matrix{
	int     m, n, l;  //< rows, cols, leading-dim
	double* vh;
	
	void Delete  ();
	void Allocate(int _m, int _n);
	void Resize  (int _m, int _n);
    void Transpose();  ///< transpose by changing layout
	
	double& operator()(int i, int j)     { assert(0 <= i && i < m && 0 <= j && j < n); return vh[l*j+i]; }
	double& operator()(int i, int j)const{ assert(0 <= i && i < m && 0 <= j && j < n); return vh[l*j+i]; }

	Matrix SubMatrix(int row, int col, int _m, int _n);
	Vector Col(int col);

	 Matrix();
	 Matrix(int _m, int _n);
	~Matrix();
};

std::ostream& operator<<(std::ostream& os, Vector& v);
std::ostream& operator<<(std::ostream& os, Matrix& m);

inline void cross3_mat(const Eigen::Vector3d& r, Eigen::Matrix3d& y){
	y(0,0) =  0.0  ; y(0,1) = -r.z(); y(0,2) =  r.y();
	y(1,0) =  r.z(); y(1,1) =  0.0  ; y(1,2) = -r.x();
	y(2,0) = -r.y(); y(2,1) =  r.x(); y(2,2) =  0.0;
}

inline void  mat_copy (const Eigen::Matrix3d& m1, Matrix&& y){
	const double* col0 = &m1(0,0);
	double*       col1 = y.vh;
	for(int j = 0; j < 3; j++, col0 += 3, col1 += y.l){
		const double* v0 = col0;
		double*       v1 = col1;
		for(int i = 0; i < 3; i++){
			*v1++ = (double)*v0++;
		}
	}
}
inline void  mat_copy (const Eigen::Matrix3d& m1, Matrix& y){
	mat_copy(m1, std::move(y));
}

template<class T, int N>
void vec_copy(const Eigen::Matrix<T, N, 1>& v, Vector&& y){
	const T* v0 = &v[0];
	double*  v1 = y.vh;
	for(int i = 0; i < N; i++)
		*v1++ = (double)*v0++;
}

template<class T, int R, int C>
void  mat_copy (const Eigen::Matrix<T, R, C>& m1, Matrix&& y){
	const T* col0 = &m1(0,0);
	double*  col1 = y .vh;
	for(int j = 0; j < C; j++, col0 += R, col1 += y.l){
		const double* v0 = col0;
		double*       v1 = col1;
		for(int i = 0; i < R; i++){
			*v1++ = (double)*v0++;
		}
	}
}

inline void vec_clear(Vector&& y, double val = 0.0){
	double* vh = y.vh;
	for(int i = 0; i < y.n; i++)
		*vh++ = val;
}
inline void vec_clear(Vector& y, double val = 0.0){
	vec_clear((Vector&&)std::move(y), val);
}
inline double vec_norm(const Vector& v){
	double vn = 0.0;
	for(int j = 0; j < v.n; j++)
		vn += (v(j)*v(j));

	return sqrt(vn);
}
inline void mat_clear(Matrix&& y, double val = 0.0){
	double* vh0 = y.vh;
	for(int j = 0; j < y.n; j++, vh0 += y.l){
		double* vh = vh0;
		for(int i = 0; i < y.m; i++, vh++){
			*vh = val;
		}
	}
}
inline void mat_clear(Matrix& y, double val = 0.0){
	mat_clear((Matrix&&)std::move(y), val);
}

inline void vec_copy(const Vector& v1, Vector&& y){
	assert(y.n == v1.n);
	memcpy(y.vh, v1.vh, sizeof(double)*y.n);
}
inline void vec_copy(const Vector& v1, Vector&& y, double k){
	assert(y.n == v1.n);
	double* vh0 = v1.vh;
	double* vh1 = y .vh;
	for(int i = 0; i < v1.n; i++)
		*vh1++ = k*(*vh0++);
}
inline void vec_copy(const Vector& v1, Vector& y){
	vec_copy(v1, (Vector&&)std::move(y));
}
inline void vec_copy(const Vector& v1, Vector& y, double k){
	vec_copy(v1, (Vector&&)std::move(y), k);
}

inline void mat_copy(const Matrix& m1, Matrix&& y){
	assert(y.m == m1.m && y.n == m1.n);
	double* col0 = m1.vh;
	double* col1 = y .vh;
	for(int j = 0; j < m1.n; j++, col0 += m1.l, col1 += y.l){
		double* v0 = col0;
		double* v1 = col1;
		memcpy(v1, v0, sizeof(double)*m1.m);
		//for(int i = 0; i < m1.m; i++){
		//	*v1++ = *v0++;
		//}
	}
}
inline void mat_copy(const Matrix& m1, Matrix&& y, double k){
	assert(y.m == m1.m && y.n == m1.n);
	double* col0 = m1.vh;
	double* col1 = y .vh;
	for(int j = 0; j < m1.n; j++, col0 += m1.l, col1 += y.l){
		double* v0 = col0;
		double* v1 = col1;
		for(int i = 0; i < m1.m; i++){
			*v1++ = k*(*v0++);
		}
	}
}
inline void mattr_copy(const Matrix& m1, Matrix&& y){
	assert(y.m == m1.n && y.n == m1.m);
	double* col0 = m1.vh;
	double* row1 = y .vh;
	for(int j = 0; j < m1.n; j++, col0 += m1.l, row1++){
		double* v0 = col0;
		double* v1 = row1;
		for(int i = 0; i < m1.m; i++){
			*v1 = *v0++;
			v1 += y.l;
		}
	}
}

inline void mat_copy(const Matrix& m1, Matrix& y){
	mat_copy(m1, (Matrix&&)std::move(y));
}
inline void mat_copy(const Matrix& m1, Matrix& y, double k){
	mat_copy(m1, (Matrix&&)std::move(y), k);
}
inline void mattr_copy(const Matrix& m1, Matrix& y){
	mattr_copy(m1, (Matrix&&)std::move(y));
}

inline void vec_add(const Vector& v1, Vector&& y){
	assert(y.n == v1.n);
	double* vh0 = v1.vh;
	double* vh1 = y .vh;
	for(int i = 0; i < v1.n; i++)
		*vh1++ += *vh0++;
}
inline void vec_add(const Vector& v1, Vector&& y, double k){
	assert(y.n == v1.n);
	double* vh0 = v1.vh;
	double* vh1 = y .vh;
	for(int i = 0; i < v1.n; i++)
		*vh1++ += k*(*vh0++);
}
inline void vec_add(const Vector& v1, Vector& y){
	vec_add(v1, (Vector&&)std::move(y));
}
inline void vec_add(const Vector& v1, Vector& y, double k){
	vec_add(v1, (Vector&&)std::move(y), k);
}

// y += m1
inline void mat_add(const Matrix& m1, Matrix&& y){
	assert(y.m == m1.m && y.n == m1.n);
	double* col0 = m1.vh;
	double* col1 = y .vh;
	for(int j = 0; j < m1.n; j++, col0 += m1.l, col1 += y.l){
		double* v0 = col0;
		double* v1 = col1;
		for(int i = 0; i < m1.m; i++){
			*v1++ += *v0++;
		}
	}
}
// y += k*m1
inline void mat_add(const Matrix& m1, Matrix&& y, double k){
	assert(y.m == m1.m && y.n == m1.n);
	double* col0 = m1.vh;
	double* col1 = y .vh;
	for(int j = 0; j < m1.n; j++, col0 += m1.l, col1 += y.l){
		double* v0 = col0;
		double* v1 = col1;
		for(int i = 0; i < m1.m; i++){
			*v1++ += k*(*v0++);
		}
	}
}
// y += m1^T
inline void mattr_add(const Matrix& m1, Matrix&& y){
	assert(y.m == m1.n && y.n == m1.m);
	double* col0 = m1.vh;
	double* row1 = y .vh;
	for(int j = 0; j < m1.n; j++, col0 += m1.l, row1++){
		double* v0 = col0;
		double* v1 = row1;
		for(int i = 0; i < m1.m; i++){
			*v1 += *v0++;
			v1 += y.l;
		}
	}
}
inline void mat_add(const Matrix& m1, Matrix& y){
	mat_add(m1, (Matrix&&)std::move(y));
}
inline void mat_add(const Matrix& m1, Matrix& y, double k){
	mat_add(m1, (Matrix&&)std::move(y), k);
}
inline void mattr_add(const Matrix& m1, Matrix& y){
	mattr_add(m1, (Matrix&&)std::move(y));
}

inline double mat_abs(const Matrix& m){
	// assumes m.m == m.l
	return cblas_dasum(m.m*m.n, m.vh, 1);
}

inline double vec_dot(const Vector& v1, const Vector& v2){
	return cblas_ddot(v1.n, v1.vh, 1, v2.vh, 1);
}

inline double quadform(const Matrix& m, const Vector& v1, const Vector& v2){
	assert(v1.n == m.m && v2.n == m.n);
	double y = 0.0;
	for(int i = 0; i < v1.n; i++)for(int j = 0; j < v2.n; j++)
		y += m(i,j)*v1(i)*v2(j);

	return y;
}

inline void mat_vec_mul(const Matrix& m1, const Vector& v, Vector&& y, double alpha, double beta){
	Vector tmp;
	if(v.vh == y.vh){
		if(tmp.n != y.n)
			tmp.Allocate(y.n);
		vec_copy(y, tmp);
		mat_vec_mul(m1, v, std::move(tmp), alpha, beta);
		vec_copy(tmp, y);
	}
	else{
		cblas_dgemv(CblasColMajor, CblasNoTrans, m1.m, m1.n, alpha, m1.vh, m1.l, v.vh, 1, beta, y.vh, 1);
	}
}
inline void mat_vec_mul(const Matrix& m1, const Vector& v, Vector& y, double alpha, double beta){
	mat_vec_mul(m1, v, (Vector&&)std::move(y), alpha, beta);
}

inline void mattr_vec_mul(const Matrix& m1, const Vector& v, Vector&& y, double alpha, double beta){
	Vector tmp;
	if(v.vh == y.vh){
		if(tmp.n != y.n)
			tmp.Allocate(y.n);
		vec_copy(y, tmp);
		mattr_vec_mul(m1, v, std::move(tmp), alpha, beta);
		vec_copy(tmp, y);
	}
	else{
		cblas_dgemv(CblasColMajor, CblasTrans, m1.m, m1.n, alpha, m1.vh, m1.l, v.vh, 1, beta, y.vh, 1);
	}
}
inline void mattr_vec_mul(const Matrix& m1, const Vector& v, Vector& y, double alpha, double beta){
	mattr_vec_mul(m1, v, (Vector&&)std::move(y), alpha, beta);
}

inline void symmat_vec_mul(const Matrix& m1, const Vector& v, Vector& y, double alpha, double beta){
	cblas_dsymv(CblasColMajor, CblasUpper, m1.m, alpha, m1.vh, m1.l, v.vh, 1, beta, y.vh, 1);
}

inline void mat_mat_mul(const Matrix& m1, const Matrix& m2, Matrix&& y, double alpha, double beta){
	Matrix tmp;
	if(m1.vh == y.vh || m2.vh == y.vh){
		if(tmp.m != y.m || tmp.n != y.n)
			tmp.Allocate(y.m, y.n);
		mat_copy(y, tmp);
		mat_mat_mul(m1, m2, std::move(tmp), alpha, beta);
		mat_copy(tmp, y);
	}
	else{
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m1.m, m2.n, m2.m, alpha, m1.vh, m1.l, m2.vh, m2.l, beta, y.vh, y.l);
	}
}
inline void mat_mat_mul(const Matrix& m1, const Matrix& m2, Matrix& y, double alpha, double beta){
	mat_mat_mul(m1, m2, (Matrix&&)std::move(y), alpha, beta);
}

inline void mattr_mat_mul(const Matrix& m1, const Matrix& m2, Matrix&& y, double alpha, double beta){
	Matrix tmp;
	if(m1.vh == y.vh || m2.vh == y.vh){
		if(tmp.m != y.m || tmp.n != y.n)
			tmp.Allocate(y.m, y.n);
		mat_copy(y, tmp);
		mattr_mat_mul(m1, m2, std::move(tmp), alpha, beta);
		mat_copy(tmp, y);
	}
	else{
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m1.n, m2.n, m1.m, alpha, m1.vh, m1.l, m2.vh, m2.l, beta, y.vh, y.l);
	}
}
inline void mattr_mat_mul(const Matrix& m1, const Matrix& m2, Matrix& y, double alpha, double beta){
	mattr_mat_mul(m1, m2, (Matrix&&)std::move(y), alpha, beta);
}

inline void mat_mattr_mul(const Matrix& m1, const Matrix& m2, Matrix& y, double alpha, double beta){
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m1.m, m2.m, m1.n, alpha, m1.vh, m1.l, m2.vh, m2.l, beta, y.vh, y.l);
}

inline void symmat_mat_mul(const Matrix& m1, const Matrix& m2, Matrix&& y, double alpha, double beta){
	cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper, m1.m, m2.n, alpha, m1.vh, m1.l, m2.vh, m2.l , beta, y.vh, y.l);
}
inline void symmat_mat_mul(const Matrix& m1, const Matrix& m2, Matrix& y, double alpha, double beta){
	symmat_mat_mul(m1, m2, (Matrix&&)std::move(y), alpha, beta);
}

inline void  mat_eye(Matrix& m){
	for(int j = 0; j < m.n; j++)for(int i = 0; i < m.m; i++)
		m(i,j) = (i == j ? 1.0f : 0.0f);
}
inline void mat_identity(Matrix& m){
	mat_eye(m);
}
inline void  mat_diag(const Vector& v, Matrix& m){
	m.Allocate(v.n, v.n);
	for(int j = 0; j < m.n; j++)for(int i = 0; i < m.m; i++)
		m(i,j) = (i == j ? v.vh[i] : 0.0f);
}

inline void  mat_eig(const Matrix& m, Vector& wr, Vector& wi, Matrix& vl, Matrix& vr){
	wr.Allocate(m.m);
	wi.Allocate(m.m);
	vl.Allocate(m.m, m.m);
	vr.Allocate(m.m, m.m);
	Matrix tmp;
	tmp.Allocate(m.m, m.m);
	mat_copy(m, tmp);
	LAPACKE_dgeev(LAPACK_COL_MAJOR, 'V', 'V', m.m, &tmp.vh[0], tmp.l, &wr.vh[0], &wi.vh[0], &vl.vh[0], vl.l, &vr.vh[0], vr.l);
}

// inverse of general matrix
inline void mat_inv_gen(const Matrix& m, Matrix& y){
	mat_copy(m, y);
	std::vector<MKL_INT> pivot; pivot.resize(m.m);
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, m.m, m.m, &y.vh[0], y.l, &pivot[0]);
	LAPACKE_dgetri(LAPACK_COL_MAJOR, m.m,      &y.vh[0], y.l, &pivot[0]);
}

// inverse of symmetric positive definite matrix
inline void mat_inv_pd(const Matrix& m, Matrix& y){
	mat_copy(m, y);
	LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', m.m, y.vh, y.l);
	LAPACKE_dpotri(LAPACK_COL_MAJOR, 'U', m.m, y.vh, y.l);
}

}
