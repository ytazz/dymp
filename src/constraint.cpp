#include <constraint.h>
#include <variable.h>
#include <solver.h>
#include <link.h>

namespace dymp{;

Constraint::Constraint(Solver* solver, int n, ID _id, int _type, real_t _scale):ID(_id){
    type       = _type;
	nelem	   = n;
	nelem_var  = 0;
	level      = 0;
	index      = 0;
	enabled    = true;
	active     = true;
	weight     = vec3_t(1.0, 1.0, 1.0);
    barrier_margin = 0.001;
	scale      = _scale;
	scale2     = scale * scale;
	scale_inv  = 1.0 / scale;
	scale2_inv = scale_inv * scale_inv;
	corrRate   = 0.1;
	corrMax    = FLT_MAX;

	solver->AddCon(this);
}

SLink* Constraint::AddSLink(Variable* var, real_t coef){
	assert(nelem == var->nelem);
	SLink* link = new SLink(var, this, coef);
	link->Connect();
	solver->links.push_back(unique_ptr<Link>(link));
	return link;
}

C2Link* Constraint::AddC2Link(Variable* var){
	assert(nelem == 2 && var->nelem == 1);
	C2Link* link = new C2Link(var, this);
	link->Connect();
	solver->links.push_back(unique_ptr<Link>(link));
	return link;
}

R2Link* Constraint::AddR2Link(Variable* var){
	assert(nelem == 1 && var->nelem == 2);
	R2Link* link = new R2Link(var, this);
	link->Connect();
	solver->links.push_back(unique_ptr<Link>(link));
	return link;
}

M2Link* Constraint::AddM2Link(Variable* var){
	assert(nelem == 2 && var->nelem == 2);
	M2Link* link = new M2Link(var, this);
	link->Connect();
	solver->links.push_back(unique_ptr<Link>(link));
	return link;
}

X3Link* Constraint::AddX3Link(Variable* var){
	assert(nelem == 3 && var->nelem == 3);
	X3Link* link = new X3Link(var, this);
	link->Connect();
	solver->links.push_back(unique_ptr<Link>(link));
	return link;
}

C3Link* Constraint::AddC3Link(Variable* var){
	assert(nelem == 3 && var->nelem == 1);
	C3Link* link = new C3Link(var, this);
	link->Connect();
	solver->links.push_back(unique_ptr<Link>(link));
	return link;
}

R3Link* Constraint::AddR3Link(Variable* var){
	assert(nelem == 1 && var->nelem == 3);
	R3Link* link = new R3Link(var, this);
	link->Connect();
	solver->links.push_back(unique_ptr<Link>(link));
	return link;
}

M3Link* Constraint::AddM3Link(Variable* var){
	assert(nelem == 3 && var->nelem == 3);
	M3Link* link = new M3Link(var, this);
	link->Connect();
	solver->links.push_back(unique_ptr<Link>(link));
	return link;
}

void Constraint::CalcError(){
	if(enabled){
		CalcDeviation();

        if( type == Type::Equality ||
            type == Type::InequalityPenalty ){
            // quadratic cost
		    for(int k = 0; k < nelem; k++)
			    e[k] = 0.5 * weight[k]*weight[k]*scale2_inv*y[k]*y[k];
        }
        if( type == Type::InequalityBarrier ){
            // logarithmic cost
			// T.B.D. take into account weight and scale
            for(int k = 0; k < nelem; k++){
			    e[k] = -solver->param.complRelaxation*log(std::max(barrier_margin, y[k]));
                e[k] = std::max(0.0, e[k]);
            }
        }
	}
	else{
		y.setZero();
		e.setZero();
	}
}

void Constraint::CalcDeviation(){
	y.setZero();
	for(Link* link : links)
		link->AddError();
}

void Constraint::RegisterCorrection(Vector&& dydvec, const vec3_t& _weight){
	for(int i = 0; i < nelem; i++)
		//dydvec[index+i] = weight[i] * dyd[i];
		dydvec(i) = (_weight[i]*scale_inv) * dyd[i];
}

void Constraint::RegisterDeviation(Vector&& yvec){
	for(int i = 0; i < nelem; i++)
		//yvec[index+i] = weight[i] * y[i];
		//yvec[offset+i] = (weight[i]*scale_inv) * y[i];
		yvec(i) = (scale_inv) * y[i];
}

void Constraint::ResetState(){
	dy.setZero();
	l .setZero();
	dl.setZero();
}

void Constraint::CalcCorrection(){
	if( type == Type::Equality ||
        type == Type::InequalityPenalty ){

        // negative gradient of quadratic cost times correction rate
	    dyd = -corrRate * y;
    }
    if( type == Type::InequalityBarrier ){
        dyd =  corrRate * y;
        //for(int k = 0; k < nelem; k++){
        //    // negative gradient of logarithmic cost times correction rate
        //    dyd[k] = corrRate * solver->param.complRelaxation/y[k];
        //}
    }

	// ただし修正幅は上限を超えないようにする
    const real_t dyd_lim = corrMax;
	real_t dyd_max = 0.0;
	for(int k = 0; k < nelem; k++)
		dyd_max = std::max(dyd_max, std::abs(dyd[k]));
	
	if(dyd_max > dyd_lim)
		dyd *= (dyd_lim / dyd_max);
    
}

//-------------------------------------------------------------------------------------------------

FixConS::FixConS(Solver* solver, ID id, SVar* var, real_t _scale):Constraint(solver, 1, id, Constraint::Type::Equality, _scale){
	desired = 0.0;
	AddSLink(var, 1.0);
}

void FixConS::CalcDeviation(){
	y[0] = ((SVar*)links[0]->var)->val - desired;
}

//-------------------------------------------------------------------------------------------------

FixConV2::FixConV2(Solver* solver, ID id, V2Var* var, real_t _scale):Constraint(solver, 2, id, Constraint::Type::Equality, _scale){
	AddSLink(var, 1.0);
}

void FixConV2::CalcDeviation(){
	y[0] = ((V2Var*)links[0]->var)->val[0] - desired[0];
	y[1] = ((V2Var*)links[0]->var)->val[1] - desired[1];
}

//-------------------------------------------------------------------------------------------------

FixConV3::FixConV3(Solver* solver, ID id, V3Var* var, real_t _scale):Constraint(solver, 3, id, Constraint::Type::Equality, _scale){
	AddSLink(var, 1.0);
}

void FixConV3::CalcDeviation(){
	y = ((V3Var*)links[0]->var)->val - desired;
}

//-------------------------------------------------------------------------------------------------

FixConQ::FixConQ(Solver* solver, ID id, QVar* var, real_t _scale):Constraint(solver, 3, id, Constraint::Type::Equality, _scale){
	AddSLink(var, 1.0);
}

void FixConQ::CalcDeviation(){
	quat_t q0 = desired;
	quat_t q1 = ((QVar*)links[0]->var)->val;
	AngleAxisd qerror(q0.conjugate() * q1);
	vec3_t axis   = qerror.axis ();
	real_t theta  = qerror.angle();
	if(theta > _pi)
		theta -= 2*_pi;
	y = q0 * (theta * axis);
}

//-------------------------------------------------------------------------------------------------

RangeConS::RangeConS(Solver* solver, ID id, SVar* var, real_t _scale):Constraint(solver, 1, id, Constraint::Type::InequalityPenalty, _scale){
	AddSLink(var, 1.0);
	real_t inf = numeric_limits<real_t>::max();
	_min = -inf;
	_max =  inf;
	on_lower = false;
	on_upper = false;
}

void RangeConS::CalcDeviation(){
	real_t s = ((SVar*)links[0]->var)->val;
	on_lower = (s <= _min);
	on_upper = (s >= _max);
	active = on_lower | on_upper;
	if(on_lower)
		y[0] = s - _min;
	if(on_upper)
		y[0] = s - _max;
}

}
