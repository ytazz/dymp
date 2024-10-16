#pragma once

#include <dymp/link.h>
#include <dymp/blas.h>
#include <dymp/id.h>

#include <map>

using namespace std;

namespace dymp{;

class Solver;
class SLink;
class C2Link;
class R2Link;
class M2Link;
class X3Link;
class C3Link;
class R3Link;
class M3Link;
class SVar;
class V2Var;
class V3Var;
class QVar;

/**
	constraint base class
 */
class Constraint : public ID{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	struct Type{
        enum{
            Equality,
            InequalityPenalty,
            InequalityBarrier,
        };
    };

	Solver*		solver;
    int         type;       ///< constraint type
	Links		links;		///< links to constrained variables
	Links       links_active;
	int         nelem;      ///< dimension of this constraint
	int         nelem_var;  ///< total dimension of variables linked with this constraint
	int         level;      ///< priority level
	int         index;      ///< index of this constraint in entire constraint variable y
	bool		enabled;	///< enabled constraint (controlled by user)
	bool		active;		///< active constraint  (for penalty inequality constraints)
	vec3_t      weight;     ///< weight of constraint error (component-wise)
	real_t      barrier_margin;   ///< margin from zero (for barrier inequality)
    real_t		scale, scale2, scale_inv, scale2_inv;		///< scaling coefficient
	
	/** error correction rate
		*/
	real_t	corrRate;

	/** maximum error correction
	 */
	real_t	corrMax;

	vec3_t		e;			///< error value
	vec3_t		y;			///< constraint error
	vec3_t		dy;			///< change of constraint error
	vec3_t		dyd;		///< desired change of constraint error
	vec3_t		l;			///< multiplier
	vec3_t		dl;			///< change of multiplier
	vec3_t		J, Jinv;	///< square sum of Jacobian row and its inverse

public:	
	SLink*		AddSLink (Variable* var, real_t coef = 1.0);
	C2Link*		AddC2Link(Variable* var);
	R2Link*		AddR2Link(Variable* var);
	M2Link*		AddM2Link(Variable* var);
	X3Link*		AddX3Link(Variable* var);
	C3Link*		AddC3Link(Variable* var);
	R3Link*		AddR3Link(Variable* var);
	M3Link*		AddM3Link(Variable* var);

	/// reset internal variables
	void ResetState();

	/// preparation
	void CalcError();
	void CalcCorrection();
	void RegisterCorrection(Vector&& dydvec, const vec3_t& _weight);
	void RegisterDeviation (Vector&& yvec);

public:
	/// virtual functions to be overridden by derived classes ///
	virtual void CalcCoef         (){}
	virtual void CalcDeviation    ();
	
public:
	
	Constraint(Solver* solver, int n, ID _id, int type, real_t _scale);
};

/** fixation constraint
 */
struct FixConS : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	real_t	desired;
	virtual void CalcDeviation();
	FixConS(Solver* solver, ID id, SVar* var, real_t _scale);
};
struct FixConV2 : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	vec2_t	desired;
	virtual void CalcDeviation();
	FixConV2(Solver* solver, ID id, V2Var* var, real_t _scale);
};
struct FixConV3 : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	vec3_t	desired;
	virtual void CalcDeviation();
	FixConV3(Solver* solver, ID id, V3Var* var, real_t _scale);
};
struct FixConQ : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	quat_t	desired;
	virtual void CalcDeviation();
	FixConQ(Solver* solver, ID id, QVar* var, real_t _scale);
};

/**	range constraint for scalar variables
 */
struct RangeConS : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	real_t	_min, _max;
	bool	on_lower, on_upper;

	virtual void CalcDeviation();
	
	RangeConS(Solver* solver, ID id, SVar* var, real_t _scale);
};

}
