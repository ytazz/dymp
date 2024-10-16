#include <dymp/variable.h>
#include <dymp/solver.h>

namespace dymp{;

Variable::Variable(int _type, Solver* solver, ID _id, real_t _scale):ID(_id){
	solver->AddVar(this);
	SetScale(_scale);

	type = _type;
	switch(type){
	case Scalar:	nelem = 1; break;
	case Vec2:      nelem = 2; break;
	case Vec3:		nelem = 3; break;
	case Quat:
	default:		nelem = 3; break;
	}

	index          = 0;
	dmax           = 1.0;
	weight         = vec3_t(0.0, 0.0, 0.0);
	index_weighted = -1;
	locked         = false;
}

void Variable::Lock(bool on){
	locked = on;
}

void Variable::SetScale(real_t sc){
	scale      = sc;
	scale2     = sc*sc;
	scale_inv  = (real_t)1.0/sc;
	scale_inv2 = scale_inv*scale_inv;
}

void Variable::ResetState(){
	dx.setZero();
	dmax2 = dmax*dmax;
}

void Variable::RegisterDelta(const Vector& dxvec){
	for(int n = 0; n < nelem; n++)
		dx[n] = dxvec(index+n);
}

}
