#pragma once

#include <dymp/types.h>
#include <dymp/blas.h>

#include <memory>

namespace dymp{;

class Variable;
class Constraint;
class Link;

typedef std::vector< Link* >	Links;

/**
	Link: Jacobian matrix from variable to constraint error
 */
class Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Variable*		var;
	Constraint*		con;
	int             index;  ///< index of var in all variables linked with con

	void Connect();

    virtual void AddError() = 0;
	virtual void RegisterCoef(Matrix&& J, const vec3_t& w) = 0;

	Link(Variable* v, Constraint* c);
};

/** link with scalar coefficient
	y = c(x) = k * x
	k is scalar
	x is scalar or vector3
	y is scalar or vector3
 */
class SLink : public Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	real_t coef, coefsqr;
	
	void SetCoef(real_t k);

	virtual void AddError    ();
	virtual void RegisterCoef(Matrix&& J, const vec3_t& w);

	SLink(Variable* v, Constraint* c, real_t k);
};

class V2Link : public Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	vec2_t	coef, coefsqr;

	void SetCoef(const vec2_t& k);

	V2Link(Variable* v, Constraint* c):Link(v, c){}
};

class V3Link : public Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	vec3_t	coef, coefsqr;

	void SetCoef(const vec3_t& k);

	V3Link(Variable* v, Constraint* c):Link(v, c){}
};

/** link with cross product matrix
	y = c(x) = k % x
	k is vector3
	x is vector3
	y is vector3
 */
class X3Link : public V3Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	virtual void AddError    ();
	virtual void RegisterCoef(Matrix&& J, const vec3_t& w);

	X3Link(Variable* v, Constraint* c):V3Link(v, c){}
};

/** link between vec2/3 constraint and scalar variable
	y = c(x) = k * x
	k is vector2/3
	x is scalar
	y is vector2/3

	* C stands for column
 */
class C2Link : public V2Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
public:
	virtual void AddError    ();
	virtual void RegisterCoef(Matrix&& J, const vec3_t& w);

	C2Link(Variable* v, Constraint* c):V2Link(v, c){}
};
class C3Link : public V3Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
public:
	virtual void AddError    ();
	virtual void RegisterCoef(Matrix&& J, const vec3_t& w);

	C3Link(Variable* v, Constraint* c):V3Link(v, c){}
};

/** link between scalar constraint and vec2/3 variable
	y = c(x) = k * x	(inner product)

	R stands for row
 */
class R2Link : public V2Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	virtual void AddError    ();
	virtual void RegisterCoef(Matrix&& J, const vec3_t& w);

	R2Link(Variable* v, Constraint* c):V2Link(v, c){}
};
class R3Link : public V3Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	virtual void AddError    ();
	virtual void RegisterCoef(Matrix&& J, const vec3_t& w);

	R3Link(Variable* v, Constraint* c):V3Link(v, c){}
};

/** general linear map
	y = A x
 */
class M2Link : public Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	mat2_t	coef;
	mat2_t	coefsqr;

public:
	void SetCoef(const mat2_t& m);

	virtual void AddError    ();
	virtual void RegisterCoef(Matrix&& J, const vec3_t& w);
	
	M2Link(Variable* v, Constraint* c):Link(v, c){}
};

class M3Link : public Link{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	mat3_t	coef;
	mat3_t	coefsqr;

public:
	void SetCoef(const mat3_t& m);

	virtual void AddError    ();
	virtual void RegisterCoef(Matrix&& J, const vec3_t& w);

	M3Link(Variable* v, Constraint* c):Link(v, c){}
};

}
