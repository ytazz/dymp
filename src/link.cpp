#include <link.h>
#include <variable.h>
#include <constraint.h>

namespace dimp3{;

Link::Link(Variable* v, Constraint* c){
	var = v;
	con = c;
	index = 0;
}

void Link::Connect(){
	var->links.push_back(this);
	con->links.push_back(this);
	index = con->nelem_var;
	con->nelem_var += var->nelem;
}

//-------------------------------------------------------------------------------------------------

SLink::SLink(Variable* v, Constraint* c, real_t k):Link(v, c){
	SetCoef(k);
}

void SLink::SetCoef(real_t k){
	coef = k;
}

void SLink::AddError(){
	if(con->nelem == 1){
		con->y[0] += coef * ((SVar *)var)->val;
	}
	else if(con->nelem == 2){
		con->y[0] += coef * ((V2Var*)var)->val[0];
		con->y[1] += coef * ((V2Var*)var)->val[1];
	}
	else{
		con->y += coef * ((V3Var*)var)->val;
	}
}

void SLink::RegisterCoef(Matrix&& J, vec3_t w){
	w *= (con->scale_inv*var->scale);
	for(int i = 0; i < con->nelem; i++)
		J(i,i) = w[i]*coef;
}	

//-------------------------------------------------------------------------------------------------

void V2Link::SetCoef(vec2_t k){
	coef = k;
}

//-------------------------------------------------------------------------------------------------

void V3Link::SetCoef(vec3_t k){
	coef = k;
}

//-------------------------------------------------------------------------------------------------

void X3Link::AddError(){
	con->y += coef.cross( ((V3Var*)var)->val );
}

void X3Link::RegisterCoef(Matrix&& J, vec3_t w){
	w *= (con->scale_inv*var->scale);
	J(0,0) =  0.0         ; J(0,1) = -w[0]*coef[2]; J(0,2) =  w[0]*coef[1];
	J(1,0) =  w[1]*coef[2]; J(1,1) =  0.0         ; J(1,2) = -w[1]*coef[0];
	J(2,0) = -w[2]*coef[1]; J(2,1) =  w[2]*coef[0]; J(2,2) =  0.0         ;
}	

//-------------------------------------------------------------------------------------------------

void C2Link::AddError(){
	con->y[0] += coef[0] * ((SVar*)var)->val;
	con->y[1] += coef[1] * ((SVar*)var)->val;
}

void C2Link::RegisterCoef(Matrix&& J, vec3_t w){
	w *= (con->scale_inv*var->scale);
	J(0,0) = w[0]*coef[0];
	J(1,0) = w[1]*coef[1];
}	

//-------------------------------------------------------------------------------------------------

void C3Link::AddError(){
	con->y += coef * dynamic_cast<SVar*>(var)->val;
}

void C3Link::RegisterCoef(Matrix&& J, vec3_t w){
	w *= (con->scale_inv*var->scale);
	J(0,0) = w[0]*coef[0];
	J(1,0) = w[1]*coef[1];
	J(2,0) = w[2]*coef[2];
}	

//-------------------------------------------------------------------------------------------------

void R2Link::AddError(){
	con->y[0] += coef.dot( dynamic_cast<V2Var*>(var)->val );
}

void R2Link::RegisterCoef(Matrix&& J, vec3_t w){
	w *= (con->scale_inv*var->scale);
	J(0,0) = w[0]*coef[0];
	J(0,1) = w[0]*coef[1];
}	

//-------------------------------------------------------------------------------------------------

void R3Link::AddError(){
	con->y[0] += coef.dot( dynamic_cast<V3Var*>(var)->val );
}

void R3Link::RegisterCoef(Matrix&& J, vec3_t w){
	w *= (con->scale_inv*var->scale);
	J(0,0) = w[0]*coef[0];
	J(0,1) = w[0]*coef[1];
	J(0,2) = w[0]*coef[2];
}	

//-------------------------------------------------------------------------------------------------

void M2Link::SetCoef(const mat2_t& m){
	coef = m;
}

void M2Link::AddError(){
	con->y[0] += coef.row(0) * dynamic_cast<V2Var*>(var)->val;
	con->y[1] += coef.row(1) * dynamic_cast<V2Var*>(var)->val;
}

void M2Link::RegisterCoef(Matrix&& J, vec3_t w){
	w *= (con->scale_inv*var->scale);
	for(int i = 0; i < 2; i++)for(int j = 0; j < 2; j++)
		J(i,j) = w[i]*coef(i,j);
}	

//-------------------------------------------------------------------------------------------------

void M3Link::SetCoef(const mat3_t& m){
	coef = m;
}

void M3Link::AddError(){
	con->y += coef * dynamic_cast<V3Var*>(var)->val;
}

void M3Link::RegisterCoef(Matrix&& J, vec3_t w){
	w *= (con->scale_inv*var->scale);
	for(int i = 0; i < 3; i++)for(int j = 0; j < 3; j++)
		J(i,j) = w[i]*coef(i,j);
}	

}
