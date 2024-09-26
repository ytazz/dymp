#pragma once

#include <link.h>
#include <id.h>
#include <util.h>

namespace dymp{;

class Solver;

class Variable : public ID{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/// variable types
	enum{
		Scalar = 1,
		Vec2   = 2,
		Vec3   = 3,
		Quat   = 4,
	};

    Solver*	 solver;
	Links	 links;			///< links to constraints
	Links    links_active;
	bool	 locked;		///< locked or not
	int	     type;			///< variable type
	int	     nelem;			///< number of elements
	int      index;
	int      index_weighted;

	real_t	scale, scale2, scale_inv, scale_inv2;	///< scaling factor, its inverse and squared inverse

	real_t   dmax, dmax2;   ///< upper limit of delta norm
	vec3_t   weight;
	vec3_t   dx;            ///< delta
	
public:
    void Lock(bool on = true);
	void SetScale(real_t sc);
    void RegisterDelta(const Vector& dxvec);

	virtual void    ResetState();

    virtual void	Reset ()            = 0;
	virtual void	Modify(real_t rate) = 0;
	
	Variable(int type, Solver* solver, ID _id, real_t _scale);
};

template<class T>
class VariableImpl : public Variable{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	T val;
	T val_tmp;

	virtual void ResetState(){
		Variable::ResetState();
		val_tmp = val;
	}

	VariableImpl(int _type, Solver* solver, ID _id, real_t _scale):Variable(_type, solver, _id, _scale){}
};

/**
	scalar variable
 */
class SVar : public VariableImpl<double>{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	virtual void Reset(){
		val = val_tmp = 0.0;
	}
	virtual void Modify(double alpha){
		val = val_tmp + (alpha*scale) * dx[0];
	}

	SVar(Solver* solver, ID _id = ID(), real_t _scale = 1.0):VariableImpl(Variable::Scalar, solver, _id, _scale){
		Reset();
	}
};

/**
	2D vector variable
 */
class V2Var : public VariableImpl<vec2_t>{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	virtual void Reset(){
		val     = vec2_t(0.0, 0.0);
		val_tmp = vec2_t(0.0, 0.0);
	}
	virtual void Modify(double alpha){
		val[0] = val_tmp[0] + (alpha*scale) * dx[0];
		val[1] = val_tmp[1] + (alpha*scale) * dx[1];
	}

	V2Var(Solver* solver, ID _id = ID(), real_t _scale = 1.0):VariableImpl(Variable::Vec2, solver, _id, _scale){
		Reset();
	}
};

/**
	3D vector variable
 */
class V3Var : public VariableImpl<vec3_t>{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	virtual void Reset(){
		val     = zero3;
		val_tmp = zero3;
	}
	virtual void Modify(double alpha){
		val = val_tmp + (alpha*scale) * dx;
	}

	V3Var(Solver* solver, ID _id = ID(), real_t _scale = 1.0):VariableImpl(Variable::Vec3, solver, _id, _scale){
		Reset();
	}
};

/**
	quaternion variable
 */
class QVar : public VariableImpl<quat_t>{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	virtual void Reset(){
		val     = unit_quat();
		val_tmp = unit_quat();
	}
	virtual void Modify(double alpha){
		double dx_norm = dx.norm();
        if(dx_norm < eps)
             val = val_tmp;
        else val = Eigen::AngleAxisd((alpha*scale)*dx_norm, dx/dx_norm) * val_tmp;
	}

	QVar(Solver* solver, ID _id = ID(), real_t _scale = 1.0):VariableImpl(Variable::Quat, solver, _id, _scale){
		Reset();
	}
};

}
