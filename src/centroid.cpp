#include <centroid.h>
#include <world.h>
#include <solver.h>
#include <util.h>
#include <rollpitchyaw.h>

namespace dymp{;

const real_t damping = 0.0;

// validify of waypoint values
inline bool is_valid(int i){
	return i != infi;
}
inline bool is_valid(real_t v){
	return v != inf;
}
inline bool is_valid(const vec2_t& v){
	return (v.x() != inf && v.y() != inf);
}
inline bool is_valid(const vec3_t& v){
	return (v.x() != inf && v.y() != inf && v.z() != inf);
}
inline bool is_valid(const quat_t& v){
	return (v.w() != inf && v.x() != inf && v.y() != inf && v.z() != inf);
}

CentroidData::End::End(){
	pos_t   = zero3;
	pos_r   = unit_quat();
	vel_t   = zero3;
	vel_r   = zero3;
	force_t = zero3;
	force_r = zero3;
	stiff   = 0.0;
	cmp     = zero2;
	force   = zero3;
	moment  = zero3;
	pos_tc  = zero3;
	iface   = -1;
	state   = Centroid::ContactState::Free;

	t_lift = t_land = 0.0;
	pos_t_lift = pos_t_land = zero3;
	pos_r_lift = pos_r_land = unit_quat();
	pos_t_lift_local = pos_t_land_local = zero3;
	pos_r_lift_local = pos_r_land_local = unit_quat();
	vel_t_ave = vel_r_ave = zero3;
		
	pos_t_weight  = one3;
	pos_r_weight  = one3;
	vel_t_weight  = one3;
	vel_r_weight  = one3;
	stiff_weight  = 1.0;
	cmp_weight    = one2;
	moment_weight = one3;
	force_weight  = one3;	
}

CentroidData::CentroidData(){
	pos_t = zero3;
	pos_r = unit_quat();
	vel_t = zero3;
	vel_r = zero3;
	L     = zero3;
	acc_t = zero3;
	acc_r = zero3;
	time     = 0.0;
	duration = 0.0;
	
	lbar = 0.0;
	pbar = rbar = etabar = zero3;
	fsum = etasum = zero3;
	p_rhs = zero3;

	pos_t_weight = one3;
	pos_r_weight = one3;
	vel_t_weight = one3;
	//vel_r_weight = one;
	L_weight     = one3;
	time_weight     = 1.0;
	duration_weight = 1.0;
}

void CentroidData::Init(Centroid* cen){
	int nend = (int)cen->ends.size();
	int ndiv = cen->param.rotationResolution;
	
	ends.resize(nend);

	I     .resize(ndiv+1);
	Iinv  .resize(ndiv+1);
	Llocal.resize(ndiv+1);
	for(int k = 0; k <= ndiv; k++){
		I   [k] = cen->param.I;
		Iinv[k] = I[k].inverse();
		Llocal[k] = zero3;
	}
}

void CentroidData::CopyVars(CentroidData& d){
	d = *this;
}

//-------------------------------------------------------------------------------------------------
// CentroidKey

CentroidKey::CentroidKey() {
	iphase = 0;
	idiv   = 0;
}

void CentroidKey::AddVar(Solver* solver) {
	cen = (Centroid*)model;
	int nend = (int)cen->ends.size();
	
	var_pos_t = new V3Var(solver, ID(VarTag::CentroidPosT, model, name + "_pos_t"), cen->scale.pt);
	var_pos_t->weight = damping*one3;
	solver->AddStateVar(var_pos_t, tick->idx);
	
	var_pos_r = new QVar(solver, ID(VarTag::CentroidPosR, model, name + "_pos_r"), cen->scale.pr);
	var_pos_r->weight = damping*one3;
	solver->AddStateVar(var_pos_r, tick->idx);

	var_vel_t = new V3Var(solver, ID(VarTag::CentroidVelT, model, name + "_vel_t"), cen->scale.vt);
	var_vel_t->weight = damping*one3;
    solver->AddStateVar(var_vel_t, tick->idx);
	
	var_L = new V3Var(solver, ID(VarTag::CentroidMomentum, model, name + "_L"), cen->scale.L);
	var_L->weight = damping*one3;
	solver->AddStateVar(var_L, tick->idx);
	
	ends.resize(nend);
	stringstream ss, ss2;
	for(int i = 0; i < nend; i++){
		ends[i].key = this;
		
		ss.str("");
		ss << name << "_end" << i;

		ends[i].var_pos_t  = new V3Var(solver, ID(VarTag::CentroidEndPos, model, ss.str() + "_pos_t"  ), cen->scale.pt);
		ends[i].var_pos_t->weight = damping*one3;
		solver->AddStateVar(ends[i].var_pos_t, tick->idx);

		ends[i].var_pos_r  = new QVar (solver, ID(VarTag::CentroidEndPos   , model, ss.str() + "_pos_r"  ), cen->scale.pr);
		ends[i].var_pos_r->weight = damping*one3;
		solver->AddStateVar(ends[i].var_pos_r , tick->idx);

		ends[i].var_vel_t  = new V3Var(solver, ID(VarTag::CentroidEndVel, model, ss.str() + "_vel_t"  ), cen->scale.vt);
		ends[i].var_vel_t->weight = damping*one3;
		solver->AddInputVar(ends[i].var_vel_t, tick->idx);
		
		ends[i].var_vel_r  = new V3Var(solver, ID(VarTag::CentroidEndVel   , model, ss.str() + "_vel_r"  ), cen->scale.vr);
		ends[i].var_vel_r->weight = damping*one3;
		solver->AddInputVar(ends[i].var_vel_r , tick->idx);
			
		if(cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			ends[i].var_stiff  = new SVar (solver, ID(VarTag::CentroidEndStiff , model, ss.str() + "_stiff"    ), cen->scale.tinv);
			ends[i].var_stiff->weight[0] = damping;
			solver->AddInputVar(ends[i].var_stiff , tick->idx);

			for(int j = 0; j < 2; j++){
				ends[i].var_cmp[j] = new SVar(solver, ID(VarTag::CentroidEndCmp, model, ss.str() + "_cmp"      ), cen->scale.pt );
				ends[i].var_cmp[j]->weight[0] = damping;
				solver->AddInputVar(ends[i].var_cmp[j], tick->idx);
			}

			ends[i].var_moment = new V3Var(solver, ID(VarTag::CentroidEndMoment, model, ss.str() + "_moment"   ), cen->scale.pt2);
			ends[i].var_moment->weight = damping*one3;
			solver->AddInputVar(ends[i].var_moment, tick->idx);
		}
		else{
			ends[i].var_force = new V3Var(solver, ID(VarTag::CentroidEndForce, model, ss.str() + "_force"   ), cen->scale.ft);
			ends[i].var_force->weight = damping*one3;
			solver->AddInputVar(ends[i].var_force, tick->idx);

			ends[i].var_moment = new V3Var(solver, ID(VarTag::CentroidEndMoment, model, ss.str() + "_moment"   ), cen->scale.fr);
			ends[i].var_moment->weight = damping*one3;
			solver->AddInputVar(ends[i].var_moment, tick->idx);
		}
	}

    var_time     = new SVar (solver, ID(VarTag::CentroidTime    , model, name + "_time"    ), cen->scale.t );
	var_duration = new SVar (solver, ID(VarTag::CentroidDuration, model, name + "_duration"), cen->scale.t );
	var_time    ->weight[0] = damping;
	var_duration->weight[0] = damping;
	solver->AddStateVar(var_time , tick->idx);
    solver->AddInputVar(var_duration, tick->idx);

	li  .resize(nend);
	li2 .resize(nend);
	pi  .resize(nend);
	ri  .resize(nend);
	fi  .resize(nend);
	etai.resize(nend);
	
	pi_cross.resize(nend);
	ri_cross.resize(nend);
	
	int ndiv = cen->param.rotationResolution;
	t     .resize(ndiv+1);
	t2    .resize(ndiv+1);
	t3    .resize(ndiv+1);
	C     .resize(ndiv+1);
	S     .resize(ndiv+1);
	C_tau .resize(ndiv+1);
	S_tau .resize(ndiv+1);
	C_lbar.resize(ndiv+1);
	S_lbar.resize(ndiv+1);

	data.v_rhs   .resize(ndiv+1);
	data.w_rhs   .resize(ndiv+1);
	data.L_rhs   .resize(ndiv+1);
	data.q_rhs   .resize(ndiv+1);
	
	R_omega .resize(ndiv+1);

	v_p     .resize(ndiv+1);
	v_v     .resize(ndiv+1);
	v_pbar  .resize(ndiv+1);
	v_rbar  .resize(ndiv+1);
	v_lbar  .resize(ndiv+1);
	v_tau   .resize(ndiv+1);

	L_v1    .resize(ndiv+1);
	L_p     .resize(ndiv+1);
	L_v     .resize(ndiv+1);
	L_tau   .resize(ndiv+1);
	L_rbar  .resize(ndiv+1);
	L_etabar.resize(ndiv+1);
	q_L1    .resize(ndiv+1);
		
	lbar_li    .resize(nend);
	pbar_li    .resize(nend);
	pbar_pi    .resize(nend);
	pbar_vi    .resize(nend);
	rbar_li    .resize(nend);
	rbar_ri    .resize(nend);
	etabar_li  .resize(nend);
	etabar_pi  .resize(nend);
	etabar_vi  .resize(nend);
	etabar_ri  .resize(nend);
	etabar_etai.resize(nend);

	p_li.resize(nend);
	p_pi.resize(nend);
	p_vi.resize(nend);
	p_ri.resize(nend);
	p_fi.resize(nend);

	v_li.resize(nend);
	v_pi.resize(nend);
	v_vi.resize(nend);
	v_ri.resize(nend);
	v_fi.resize(nend);

	L_li.resize  (nend);
	L_pi.resize  (nend);
	L_vi.resize  (nend);
	L_ri.resize  (nend);
	L_fi.resize  (nend);
	L_etai.resize(nend);

	q_li.resize  (nend);
	q_pi.resize  (nend);
	q_vi.resize  (nend);
	q_ri.resize  (nend);
	q_fi.resize  (nend);
	q_etai.resize(nend);

	for(int iend = 0; iend < nend; iend++){
		v_li[iend].resize(ndiv+1);
		v_pi[iend].resize(ndiv+1);
		v_vi[iend].resize(ndiv+1);
		v_ri[iend].resize(ndiv+1);
		v_fi[iend].resize(ndiv+1);

		L_li  [iend].resize(ndiv+1);
		L_pi  [iend].resize(ndiv+1);
		L_vi  [iend].resize(ndiv+1);
		L_ri  [iend].resize(ndiv+1);
		L_fi  [iend].resize(ndiv+1);
		L_etai[iend].resize(ndiv+1);
	}		

	data.Init(cen);
	data_des.Init(cen);
}

void CentroidKey::AddCon(Solver* solver) {
	CentroidKey* nextObj = (CentroidKey*)next;

    int nend  = (int)cen->ends .size();
    int nface = (int)cen->faces.size();

    if(next){
		con_pos_t = new CentroidPosConT(solver, name + "_pos_t", this, cen->scale.pt);
		solver->AddTransitionCon       (con_pos_t, tick->idx);
		
		con_pos_r = new CentroidPosConR(solver, name + "_pos_r", this, cen->scale.pr);
		solver->AddTransitionCon       (con_pos_r, tick->idx);
		
		con_vel_t = new CentroidVelConT(solver, name + "_vel_t", this, cen->scale.vt);
		solver->AddTransitionCon       (con_vel_t, tick->idx);
        
		con_L = new CentroidLCon(solver, name + "_L", this, cen->scale.L);
		solver->AddTransitionCon       (con_L, tick->idx);
		
		con_time  = new CentroidTimeCon(solver, name + "_time" , this, cen->scale.t );
		solver->AddTransitionCon(con_time , tick->idx);		
        
        con_duration_range[0] = new CentroidDurationRangeCon(solver, name + "_duration", this,  1.0, cen->scale.t);
        con_duration_range[1] = new CentroidDurationRangeCon(solver, name + "_duration", this, -1.0, cen->scale.t);
        solver->AddCostCon(con_duration_range[0], tick->idx);
        solver->AddCostCon(con_duration_range[1], tick->idx);
    }

    con_des_pos_t = new FixConV3(solver, ID(ConTag::CentroidDesPosT, model, name + "_des_pos_t"), var_pos_t, cen->scale.pt);
	solver->AddCostCon(con_des_pos_t, tick->idx);

	con_des_pos_r = new FixConQ (solver, ID(ConTag::CentroidDesPosR, model, name + "_des_pos_r"), var_pos_r, cen->scale.pr);
	solver->AddCostCon(con_des_pos_r, tick->idx);
		
	con_des_vel_t = new FixConV3(solver, ID(ConTag::CentroidDesVelT, model, name + "_des_vel_t"), var_vel_t, cen->scale.vt);
    solver->AddCostCon(con_des_vel_t, tick->idx);

	con_des_L = new FixConV3(solver, ID(ConTag::CentroidDesMomentum, model, name + "_des_L"), var_L, cen->scale.L);
	solver->AddCostCon(con_des_L, tick->idx);

	con_des_time     = new FixConS (solver, ID(ConTag::CentroidDesTime, model, name + "_des_time"    ), var_time    , cen->scale.t );
	con_des_duration = new FixConS (solver, ID(ConTag::CentroidDesDuration, model, name + "_des_duration"), var_duration, cen->scale.t );
    solver->AddCostCon(con_des_time    , tick->idx);
    solver->AddCostCon(con_des_duration, tick->idx);

	stringstream ss;
	for(int i = 0; i < nend; i++){
		ss.str("");
		ss << name << "_end" << i;

        if(next){
            ends[i].con_pos_t = new CentroidEndPosConT(solver, ss.str() + "_pos_t", this, i, cen->scale.pt);
            solver->AddTransitionCon(ends[i].con_pos_t, tick->idx);

			ends[i].con_pos_r = new CentroidEndPosConR(solver, ss.str() + "_pos_r", this, i, cen->scale.pr);
			solver->AddTransitionCon(ends[i].con_pos_r, tick->idx);

			if(cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
				ends[i].con_stiff_range     = new RangeConS (solver, ID(ConTag::CentroidEndStiffRange, model, ss.str() + "_stiff_range" ), ends[i].var_stiff, cen->scale.tinv);
				solver->AddCostCon(ends[i].con_stiff_range    , tick->idx);
			}
        }
	
		for(int j = 0; j < 3; j++){
			ends[i].con_pos_range[j][0] = new CentroidEndPosRangeCon(solver, name + "_pos_range", this, i, j,  1.0, cen->scale.pt);
			ends[i].con_pos_range[j][1] = new CentroidEndPosRangeCon(solver, name + "_pos_range", this, i, j, -1.0, cen->scale.pt);
			solver->AddCostCon(ends[i].con_pos_range[j][0], tick->idx);
			solver->AddCostCon(ends[i].con_pos_range[j][1], tick->idx);
		}

		ends[i].con_des_pos_t = new FixConV3(solver, ID(ConTag::CentroidDesEndPosT, model, name + "_des_pos"), ends[i].var_pos_t, cen->scale.pt);
		ends[i].con_des_vel_t = new FixConV3(solver, ID(ConTag::CentroidDesEndVelT, model, name + "_des_vel"), ends[i].var_vel_t, cen->scale.vt);
		solver->AddCostCon(ends[i].con_des_pos_t , tick->idx);
		solver->AddCostCon(ends[i].con_des_vel_t , tick->idx);
		
		ends[i].con_des_pos_r  = new FixConQ (solver, ID(ConTag::CentroidDesEndPosR, model, name + "_des_pos"   ), ends[i].var_pos_r , cen->scale.pr  );
		ends[i].con_des_vel_r  = new FixConV3(solver, ID(ConTag::CentroidDesEndVelR, model, name + "_des_vel"   ), ends[i].var_vel_r , cen->scale.vr  );
		solver->AddCostCon(ends[i].con_des_pos_r , tick->idx);
		solver->AddCostCon(ends[i].con_des_vel_r , tick->idx);

		if(cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){	
			ends[i].con_des_stiff  = new FixConS (solver, ID(ConTag::CentroidDesEndStiff, model, name + "_des_stiff" ), ends[i].var_stiff , cen->scale.tinv);
			solver->AddCostCon(ends[i].con_des_stiff , tick->idx);

			for(int j = 0; j < 2; j++){
				ends[i].con_des_cmp[j] = new FixConS(solver, ID(ConTag::CentroidDesEndCmp, model, name + "_des_cmp"   ), ends[i].var_cmp[j]   , cen->scale.pt  );
				solver->AddCostCon(ends[i].con_des_cmp[j], tick->idx);
			}

			ends[i].con_des_moment = new FixConV3(solver, ID(ConTag::CentroidDesEndMoment, model, name + "_des_moment"), ends[i].var_moment, cen->scale.pt2 );
			solver->AddCostCon(ends[i].con_des_moment, tick->idx);
		}
		else{
			ends[i].con_des_force = new FixConV3(solver, ID(ConTag::CentroidDesEndForce, model, name + "_des_force"), ends[i].var_force, cen->scale.ft);
			solver->AddCostCon(ends[i].con_des_force, tick->idx);

			ends[i].con_des_moment = new FixConV3(solver, ID(ConTag::CentroidDesEndMoment, model, name + "_des_moment"), ends[i].var_moment, cen->scale.fr);
			solver->AddCostCon(ends[i].con_des_moment, tick->idx);
		}
		
		/// contact force constraints
		if(next){
			real_t sc;
			sc = (cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness ? cen->scale.pt : cen->scale.ft);
			ends[i].con_friction[0][0] = new CentroidEndFrictionCon(solver, name + "_friction", this, i, 0, 0, sc);
			ends[i].con_friction[0][1] = new CentroidEndFrictionCon(solver, name + "_friction", this, i, 0, 1, sc);
			ends[i].con_friction[1][0] = new CentroidEndFrictionCon(solver, name + "_friction", this, i, 1, 0, sc);
			ends[i].con_friction[1][1] = new CentroidEndFrictionCon(solver, name + "_friction", this, i, 1, 1, sc);
			solver->AddCostCon(ends[i].con_friction[0][0], tick->idx);
			solver->AddCostCon(ends[i].con_friction[0][1], tick->idx);
			solver->AddCostCon(ends[i].con_friction[1][0], tick->idx);
			solver->AddCostCon(ends[i].con_friction[1][1], tick->idx);
			
			sc = (cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness ? cen->scale.pt2 : cen->scale.fr);
			ends[i].con_moment[0][0] = new CentroidEndMomentCon  (solver, name + "_moment00", this, i, 0, 0, sc);
			ends[i].con_moment[0][1] = new CentroidEndMomentCon  (solver, name + "_moment01", this, i, 0, 1, sc);
			ends[i].con_moment[1][0] = new CentroidEndMomentCon  (solver, name + "_moment10", this, i, 1, 0, sc);
			ends[i].con_moment[1][1] = new CentroidEndMomentCon  (solver, name + "_moment11", this, i, 1, 1, sc);
			ends[i].con_moment[2][0] = new CentroidEndMomentCon  (solver, name + "_moment20", this, i, 2, 0, sc);
			ends[i].con_moment[2][1] = new CentroidEndMomentCon  (solver, name + "_moment21", this, i, 2, 1, sc);
			solver->AddCostCon(ends[i].con_moment[0][0], tick->idx);
			solver->AddCostCon(ends[i].con_moment[0][1], tick->idx);
			solver->AddCostCon(ends[i].con_moment[1][0], tick->idx);
			solver->AddCostCon(ends[i].con_moment[1][1], tick->idx);
			solver->AddCostCon(ends[i].con_moment[2][0], tick->idx);
			solver->AddCostCon(ends[i].con_moment[2][1], tick->idx);
		}
        
        ends[i].con_contact.resize(nface);
		for(int j = 0; j < nface; j++){
		    ends[i].con_contact[j] = new CentroidEndContactCon (solver, name + "_contact" , this, i, j, cen->scale.pt);
            
			solver->AddCostCon(ends[i].con_contact[j], tick->idx);
        }
	}
}

void CentroidKey::Prepare() {
	tick->time = var_time->val;

	int nend = (int)ends.size();
	
	// copy variables to data
	data.pos_t    = var_pos_t->val;
	data.pos_r    = var_pos_r->val;
	data.vel_t    = var_vel_t->val;
	data.L        = var_L    ->val;
	data.time     = var_time ->val;
	data.duration = var_duration->val;

	for(int i = 0; i < nend; i++){
		End& end = ends[i];
		CentroidData::End & dend     = data    .ends [i];
		CentroidData::End & dend_des = data_des.ends [i];
	
		dend.pos_t = end.var_pos_t->val;
		dend.pos_r = end.var_pos_r->val;
		dend.vel_t = end.var_vel_t->val;
		dend.vel_r = end.var_vel_r->val;
	
		dend.state = dend_des.state;
        dend.iface = dend_des.iface;
		if(next){
			if(cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){	
				dend.stiff  = end.var_stiff ->val;
				dend.cmp[0] = end.var_cmp[0]->val;
				dend.cmp[1] = end.var_cmp[1]->val;
			}
			else{
				dend.force = end.var_force->val;
			}
			dend.moment = end.var_moment->val;
		}
	}
	
	if(!next)
		return;

    cen = (Centroid*)model;
	const real_t eps2 = 1.0e-06;

	int ndiv = cen->param.rotationResolution;
	
	dtau = data.duration/(real_t)ndiv;
	for(int k = 0; k <= ndiv; k++){
		t [k] = k*dtau;
		t2[k] = t[k]*t[k];
		t3[k] = t2[k]*t[k];
	}

	g = vec3_t(0.0, 0.0, cen->param.g);
	m = cen->param.m;

	for(int i = 0; i < nend; i++){
		pi [i] = data.ends[i].pos_t + data.ends[i].vel_t*(data.duration/2.0);
		pi_cross[i] = cross_mat(pi[i]);

		fi  [i] = data.ends[i].force;
		etai[i] = data.ends[i].moment;
	}

	if(cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
		// to avoid singularity in flight phase
		l2sum = eps2*eps2;
		
		psum = zero3;
		rsum = zero3;
		data.etabar = zero3;
		for(int i = 0; i < nend; i++){
			li [i] = data.ends[i].stiff;
			li2[i] = li[i]*li[i];
		
			ri      [i] = vec3_t(data.ends[i].cmp[0], data.ends[i].cmp[1], 0.0);
			ri_cross[i] = cross_mat(ri[i]);
		
			l2sum   += li2[i];
			psum    += li2[i]*pi[i];
			rsum    += li2[i]*ri[i];
			data.etabar  += li2[i]*(etai[i] - pi[i].cross(ri[i]));
		}

		data.lbar = sqrt(l2sum);
		data.pbar = (psum + g)/l2sum;
		data.rbar = rsum/l2sum;
		data.etabar += l2sum*data.pbar.cross(data.rbar);

		pbar_cross = cross_mat(data.pbar);
		rbar_cross = cross_mat(data.rbar);

		for(int k = 0; k <= ndiv; k++){
			C[k]      = cosh(data.lbar*t[k]);
			S[k]      = sinh(data.lbar*t[k]);
			C_tau [k] = ((real_t)k/(real_t)ndiv)*data.lbar*S[k];
			S_tau [k] = ((real_t)k/(real_t)ndiv)*data.lbar*C[k];
			C_lbar[k] = t[k]*S[k];
			S_lbar[k] = t[k]*C[k];
		}
	
		data.p_rhs = data.pbar + data.rbar + C[ndiv]*(data.pos_t - (data.pbar + data.rbar)) + (S[ndiv]/data.lbar)*data.vel_t;

		for(int k = 0; k <= ndiv; k++){
			data.v_rhs[k] = (data.lbar*S[k])*(data.pos_t - (data.pbar + data.rbar)) + C[k]*data.vel_t;
			data.L_rhs[k] = data.L + m*( (data.v_rhs[k] - data.vel_t).cross(data.rbar) + t[k]*data.etabar );
		}
	}
	else{
		data.fsum   = zero3;
		data.etasum = zero3;
		for(int i = 0; i < nend; i++){
			data.fsum   += fi[i];
			data.etasum += pi[i].cross(fi[i]) + etai[i];
		}
		
		data.p_rhs       = data.pos_t + data.vel_t*t[ndiv] + (data.fsum/m - g)*(t[ndiv]*t[ndiv]/2.0);
		data.v_rhs[ndiv] = data.vel_t + (data.fsum/m - g)*t[ndiv];

		for(int k = 0; k <= ndiv; k++){
			data.L_rhs[k] = data.L + data.etasum*t[k] + data.fsum.cross(data.pos_t*t[k] + data.vel_t*(t2[k]/2.0) - (t3[k]/6.0)*g);
		}		
	}

	for(int k = 0; k <= ndiv; k++){
		data.w_rhs[k] = data.Iinv[k]*(data.L_rhs[k] - data.Llocal[k]);
	}

	// q[k+1] = q(w(tk+dtau*(ndiv-1)))*q(w(tk+dtau*(ndiv-2))) * ... * q(w(tk)) * q[k]
	data.q_rhs[0] = data.pos_r;
	for(int k = 0; k < ndiv; k++){
		data.q_rhs[k+1] = rot_quat(data.w_rhs[k]*dtau)*data.q_rhs[k];
	}
}

void CentroidKey::Prepare2(){
	int nend = (int)ends.size();
	for(int i = 0; i < nend; i++){
		CentroidKey* k0 = this;
		CentroidKey* k1 = this;

		// find lift-off and landing phase
		while(k0->prev && ((CentroidKey*)k0->prev)->data_des.ends[i].iface == -1)
			k0 = (CentroidKey*)k0->prev;
		while(k1->next && k1->data_des.ends[i].iface == -1)
			k1 = (CentroidKey*)k1->next;

		CentroidData& d0 = k0->data;
		CentroidData& d1 = k1->data;

		data.ends[i].t_lift = d0.time;
		data.ends[i].t_land = d1.time;
		data.ends[i].pos_t_lift       = d0.ends[i].pos_t;
		data.ends[i].pos_r_lift       = d0.ends[i].pos_r;
		data.ends[i].pos_t_land       = d1.ends[i].pos_t;
		data.ends[i].pos_r_land       = d1.ends[i].pos_r;
		data.ends[i].pos_t_lift_local = d0.pos_r.conjugate()*(d0.ends[i].pos_t - d0.pos_t);
		data.ends[i].pos_r_lift_local = d0.pos_r.conjugate()*d0.ends[i].pos_r;
		data.ends[i].pos_t_land_local = d1.pos_r.conjugate()*(d1.ends[i].pos_t - d1.pos_t);
		data.ends[i].pos_r_land_local = d1.pos_r.conjugate()*d1.ends[i].pos_r;
	}
    for(int i = 0; i < nend; i++){
		CentroidKey* km1 = (prev && ((CentroidKey*)prev)->data_des.ends[i].state == Centroid::ContactState::Free) ? (CentroidKey*)prev : this;
		CentroidKey* k1  = (next && data_des.ends[i].state == Centroid::ContactState::Free) ? (CentroidKey*)next : this;

		CentroidData& dm1 = km1->data;
		CentroidData& d1  = k1 ->data;

		vec3_t anglem1 = ToRollPitchYaw(dm1.ends[i].pos_r);
		vec3_t angle1  = ToRollPitchYaw(d1 .ends[i].pos_r);
		data.ends[i].vel_t_ave = (d1.ends[i].pos_t - dm1.ends[i].pos_t)/(d1.time - dm1.time);
		data.ends[i].vel_r_ave = (angle1 - anglem1)/(d1.time - dm1.time);
	}
}

void CentroidKey::PrepareStep(){
	int nend = (int)ends.size();
	int ndiv = cen->param.rotationResolution;

	R_omega[ndiv] = mat3_t::Identity();
	for(int k = ndiv-1; k >= 0; k--){
		mat3_t R = rot_quat(data.w_rhs[k]*dtau).toRotationMatrix();
		R_omega[k] = R_omega[k+1]*R;
	}

	q_q = R_omega[0];
	q_tau = zero3;
	for(int k = 0; k < ndiv; k++){
		mat3_t tmp = R_omega[k+1]*rot_jacobian(data.w_rhs[k]*dtau);
		q_L1[k] = tmp*data.Iinv[k]*dtau;
		q_tau  += tmp*data.w_rhs[k]/(real_t)ndiv;
	}

	// calc direct coefficients
	if(cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
		// lbar coef
		for(int i = 0; i < nend; i++){
			lbar_li[i] = li[i]/data.lbar;
		}
		// pbar coef
		pbar_lbar = -2.0*(psum + g)/(l2sum*data.lbar);
		for(int i = 0; i < nend; i++){
			pbar_li[i] = 2.0*li[i]*pi[i]/l2sum;
			pbar_pi[i] = li2[i]/l2sum;
			pbar_vi[i] = pbar_pi[i]*(data.duration/2.0);
		}
		// rbar coef
		rbar_lbar = -2.0*rsum/(l2sum*data.lbar);
		for(int i = 0; i < nend; i++){
			rbar_li[i] = 2.0*li[i]*ri[i]/l2sum;
			rbar_ri[i] = li2[i]/l2sum;
		}
		// etabar coef
		etabar_lbar =  2.0*data.lbar*data.pbar.cross(data.rbar);
		etabar_pbar = -l2sum*rbar_cross;
		etabar_rbar =  l2sum*pbar_cross;
		for(int i = 0; i < nend; i++){
			etabar_li  [i] =  2.0*li[i]*(etai[i] - pi[i].cross(ri[i]));
			etabar_pi  [i] =  li2[i]*ri_cross[i];
			etabar_vi  [i] =  etabar_pi[i]*(data.duration/2.0);
			etabar_ri  [i] = -li2[i]*pi_cross[i];
			etabar_etai[i] =  li2[i];
		}
		// p' coef
		p_C    =  data.pos_t - (data.pbar + data.rbar);
		p_S    =  (1.0/data.lbar)*data.vel_t;
		p_p    =  C[ndiv];
		p_v    =  (S[ndiv])/data.lbar;
		p_lbar = -(S[ndiv]/l2sum)*data.vel_t;
		p_pbar =  1.0 - C[ndiv];
		p_rbar =  1.0 - C[ndiv];
		// v' coef
		v_C = data.vel_t;
		v_S = data.lbar*(data.pos_t - (data.pbar + data.rbar));
		for(int k = 0; k <= ndiv; k++){
			v_p   [k] =  data.lbar*S[k];
			v_v   [k] =  C[k];
			v_lbar[k] =  S[k]*(data.pos_t - (data.pbar + data.rbar));
			v_pbar[k] = -data.lbar*S[k];
			v_rbar[k] = -data.lbar*S[k];
		}
		
		L_L    =  1;
		for(int k = 0; k <= ndiv; k++){
			L_v1    [k] = (k == 0 ? mat3_t::Zero() : mat3_t(-m*rbar_cross));
			L_p     [k] = mat3_t::Zero();
			L_v     [k] = (k == 0 ? mat3_t::Zero() : mat3_t( m*rbar_cross));
			L_rbar  [k] = m*(cross_mat(data.v_rhs[k] - data.vel_t));
			L_etabar[k] = m*t[k];
			L_tau   [k] = m*((real_t)k/(real_t)ndiv)*data.etabar;
		}

		// calc dependent coefficients

		// p' coef
		p_tau = p_C*C_tau[ndiv] + p_S*S_tau[ndiv];
		for(int i = 0; i < nend; i++){
			p_li[i] = 
				(p_C*C_lbar[ndiv] + p_S*S_lbar[ndiv] + p_lbar)*lbar_li[i] +
				 p_pbar*(pbar_lbar*lbar_li[i] + pbar_li[i]) +
				 p_rbar*(rbar_lbar*lbar_li[i] + rbar_li[i]);
			p_pi[i] = p_pbar*pbar_pi[i];
			p_vi[i] = p_pbar*pbar_vi[i];
			p_ri[i] = p_rbar*rbar_ri[i];
		}
		// v' coef
		for(int k = 0; k <= ndiv; k++){
			v_tau[k]  = v_C*C_tau[k] + v_S*S_tau[k];
		}
		for(int i = 0; i < nend; i++){
			for(int k = 0; k <= ndiv; k++){
				v_li[i][k] = 
					(v_C*C_lbar[k] + v_S*S_lbar[k] + v_lbar[k])*lbar_li[i] +
					 v_pbar[k]*(pbar_lbar*lbar_li[i] + pbar_li[i]) +
					 v_rbar[k]*(rbar_lbar*lbar_li[i] + rbar_li[i]);
				v_pi[i][k] = v_pbar[k]*pbar_pi[i];
				v_vi[i][k] = v_pbar[k]*pbar_vi[i];
				v_ri[i][k] = v_rbar[k]*rbar_ri[i];
			}
		}
		// L' coef
		for(int k = 0; k <= ndiv; k++){
			L_p  [k] += L_v1[k]*v_p  [k];
			L_v  [k] += L_v1[k]*v_v  [k];
			L_tau[k] += L_v1[k]*v_tau[k];
		}
		for(int i = 0; i < nend; i++){
			for(int k = 0; k <= ndiv; k++){
				L_li[i][k] = 
					L_v1[k]*v_li[i][k] +
					L_rbar  [k]*(rbar_li[i] + rbar_lbar*lbar_li[i]);
					L_etabar[k]*(
						etabar_li[i] + 
						etabar_lbar*lbar_li[i] +
						etabar_pbar*(pbar_lbar*lbar_li[i] + pbar_li[i]) +
						etabar_rbar*(rbar_lbar*lbar_li[i] + rbar_li[i])
						);
				L_pi[i][k] = 
					L_v1[k]*v_pi[i][k] +
					L_etabar[k]*(etabar_pi[i] + etabar_pbar*pbar_pi[i]);
				L_vi[i][k] = 
					L_v1[k]*v_vi[i][k] +
					L_etabar[k]*(etabar_vi[i] + etabar_pbar*pbar_vi[i]);
				L_ri[i][k] = 
					L_v1[k]*v_ri[i][k] +
					L_rbar  [k]*rbar_ri[i] +
					L_etabar[k]*(etabar_ri[i] + etabar_rbar*rbar_ri[i]);
				L_etai[i][k] = 
					L_etabar[k]*etabar_etai[i];
			}
		}
		// q' coef
		q_L = mat3_t::Zero();
		q_p = mat3_t::Zero();
		q_v = mat3_t::Zero();
		for(int k = 0; k < ndiv; k++){
			q_L   += q_L1[k]*L_L;
			q_p   += q_L1[k]*L_p  [k];
			q_v   += q_L1[k]*L_v  [k];
			q_tau += q_L1[k]*L_tau[k];
		}
		for(int i = 0; i < nend; i++){
			q_li  [i] = vec3_t::Zero();
			q_pi  [i] = mat3_t::Zero();
			q_vi  [i] = mat3_t::Zero();
			q_ri  [i] = mat3_t::Zero();
			q_etai[i] = mat3_t::Zero();
			for(int k = 0; k < ndiv; k++){
				q_li  [i] += q_L1[k]*L_li  [i][k];
				q_pi  [i] += q_L1[k]*L_pi  [i][k];
				q_vi  [i] += q_L1[k]*L_vi  [i][k];
				q_ri  [i] += q_L1[k]*L_ri  [i][k];
				q_etai[i] += q_L1[k]*L_etai[i][k];
			}
		}
	}
	// direct parametrization
	else{
		// p' coef
		p_p   = 1.0;
		p_v   = t[ndiv];
		p_tau = data.vel_t + (data.fsum/m - g)*t[ndiv];
		for(int i = 0; i < nend; i++){
			p_fi[i] = (1.0/m)*(t2[ndiv]/2.0);
		}
		// v' coef
		v_v[ndiv]   = 1.0;
		v_tau[ndiv] = (data.fsum/m - g);
		for(int i = 0; i < nend; i++){
			v_fi[i][ndiv] = (1.0/m)*t[ndiv];
		}
		// L' coef
		L_L = 1.0;
		for(int k = 0; k <= ndiv; k++){
			L_p  [k] = cross_mat(data.fsum)*t[k];
			L_v  [k] = cross_mat(data.fsum)*(t2[k]/2.0);
			L_tau[k] = (data.etasum + data.fsum.cross(data.pos_t + data.vel_t*t[k] - g*(t2[k]/2.0)))*((real_t)k/(real_t)ndiv);
		}
		for(int i = 0; i < nend; i++){
			for(int k = 0; k <= ndiv; k++){
				L_pi[i][k] = -cross_mat(fi[i])*t[k];
				L_vi[i][k] = -cross_mat(fi[i])*(t2[k]/2.0);
				L_fi[i][k] =  cross_mat(pi[i]*t[k] - (data.pos_t*t[k] + data.vel_t*(t2[k]/2.0) - g*(t3[k]/6.0)));
				L_etai[i][k] = t[k];
			}
		}
		// q' coef
		q_L = mat3_t::Zero();
		q_p = mat3_t::Zero();
		q_v = mat3_t::Zero();
		for(int k = 0; k < ndiv; k++){
			q_L   += q_L1[k]*L_L;
			q_p   += q_L1[k]*L_p  [k];
			q_v   += q_L1[k]*L_v  [k];
			q_tau += q_L1[k]*L_tau[k];
		}
		for(int i = 0; i < nend; i++){
			q_pi  [i] = mat3_t::Zero();
			q_vi  [i] = mat3_t::Zero();
			q_fi  [i] = mat3_t::Zero();
			q_etai[i] = mat3_t::Zero();
			for(int k = 0; k < ndiv; k++){
				q_pi  [i] += q_L1[k]*L_pi  [i][k];
				q_vi  [i] += q_L1[k]*L_vi  [i][k];
				q_fi  [i] += q_L1[k]*L_fi  [i][k];
				q_etai[i] += q_L1[k]*L_etai[i][k];
			}
		}
	}

}

void CentroidKey::Finish(){
	tick->time = var_time->val;
	
	//int ndiv = (!cen->phases.empty() ? cen->phases[iphase].ndiv : 1);
	//var_duration->val = std::min(std::max(cen->param.durationMin/ndiv, var_duration->val), cen->param.durationMax/ndiv);
	var_duration->val = std::min(std::max(cen->param.durationMin, var_duration->val), cen->param.durationMax);
	//DSTR << "time: " << var_time->val << " duration: " << var_duration->val << endl;

	real_t dmax = 0.0;
	for(int i = 0; i < ends.size(); i++){
		End& end = ends[i];
		if(cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			end.var_stiff->val = std::min(std::max(0.0, end.var_stiff->val), 100.0);
		}
	}
	
	Prepare();	
}

void CentroidKey::Draw(render::Canvas* canvas, render::Config* conf) {
	canvas->SetPointSize(5.0f);
	canvas->Point(data.pos_t);
	
    int nend = (int)ends.size();
	for(int i = 0; i < nend; i++){
		if(data_des.ends[i].iface == -1)
			continue;

		canvas->SetPointSize(2.0f);
		canvas->SetPointColor("black");
		canvas->Point(data.ends[i].pos_t);

		canvas->SetPointColor("blue");
		canvas->Point(data.ends[i].pos_t + vec3_t(data.ends[i].cmp[0], data.ends[i].cmp[1], 0.0));
	}
}

//-------------------------------------------------------------------------------------------------
// Centroid

Centroid::Param::Param() {
	g  = 9.8;
	m  = 1.0;
	I  = mat3_t::Identity();
	mu = 1.0;

    durationMin = 0.1;
	durationMax = 1.0;

    bodyRangeMin = vec3_t(-0.1, -0.1, -0.1);
    bodyRangeMax = vec3_t( 0.1,  0.1,  0.1);

    swingHeight = 0.05;
	swingSlope  = 0.0;
    swingLiftMargin = 0.0;
	swingLandMargin = 0.0;
	swingTurnMargin = 0.0;

	contactMargin = 0.0;

	complWeight = 100.0;

	enableRotation   = true;
	rotationResolution = 1;
	endInterpolation         = EndInterpolation::CycloidGlobal;
	endWrenchParametrization = EndWrenchParametrization::Stiffness;
}

//-------------------------------------------------------------------------------------------------

Centroid::End::End(){
    stiffnessMax    = 1.0;
	lockOri         = false;
	lockCmp         = false;
	lockMoment      = false;
}

//-------------------------------------------------------------------------------------------------

Centroid::Face::Face(){
	normal       = vec3_t(0.0, 0.0, 1.0);
}

//-------------------------------------------------------------------------------------------------

Centroid::Waypoint::End::Value::Value(){
	pos_t = inf3;
	pos_r = inf3;
	vel_t = inf3;
	vel_r = inf3;
	iface = infi;
}
Centroid::Waypoint::End::Value::Value(const vec3_t& _p, const vec3_t& _q, const vec3_t& _v, const vec3_t& _w, int _iface){
	pos_t = _p;
	pos_r = _q;
	vel_t = _v;
	vel_r = _w;
	iface = _iface;
}
Centroid::Waypoint::End::Weight::Weight(){
	pos_t  = inf3;
	pos_r  = inf3;
	vel_t  = inf3;
	vel_r  = inf3;
	stiff  = inf;
	cmp    = inf2;
	moment = inf3;
}
Centroid::Waypoint::End::Weight::Weight(const vec3_t& _p, const vec3_t& _q, const vec3_t& _v, const vec3_t& _w, real_t _l, const vec2_t& _r, const vec3_t& _eta){
	pos_t  = _p;
	pos_r  = _q;
	vel_t  = _v;
	vel_r  = _w;
	stiff  = _l;
	cmp    = _r;
	moment = _eta;
}
Centroid::Waypoint::End::End(){
}
Centroid::Waypoint::Value::Value(){
	time     = inf;
	duration = inf;
	pos_t    = vec3_t(inf, inf, inf);
	pos_r    = vec3_t(inf, inf, inf);
	vel_t    = vec3_t(inf, inf, inf);
	vel_r    = vec3_t(inf, inf, inf);
}
Centroid::Waypoint::Value::Value(real_t _t, real_t _tau, const vec3_t& _p, const vec3_t& _q, const vec3_t& _v, const vec3_t& _w){
	time     = _t;
	duration = _tau;
	pos_t    = _p;
	pos_r    = _q;
	vel_t    = _v;
	vel_r    = _w;
}
Centroid::Waypoint::Weight::Weight(){
	time     = inf;
	duration = inf;
	pos_t    = inf3;
	pos_r    = inf3;
	vel_t    = inf3;
	L        = inf3;
}
Centroid::Waypoint::Weight::Weight(real_t _t, real_t _tau, const vec3_t& _p, const vec3_t& _q, const vec3_t& _v, const vec3_t& _L){
	time     = _t;
	duration = _tau;
	pos_t    = _p;
	pos_r    = _q;
	vel_t    = _v;
	L        = _L;
}

Centroid::Waypoint::Waypoint() {

}

//-------------------------------------------------------------------------------------------------

Centroid::Centroid(World* w, string n) :Model(w, n) {
	callback = 0;
}

Centroid::~Centroid() {
}

void Centroid::SetScaling(){
	// moment of inertial of solid sphere = 0.4 m r^2
	// r = sqrt(I/(0.4m))
	scale.l    = 1.0;//sqrt(10*param.I[0][0]/(0.4*param.m));//0.5;  //< unit length
	scale.t    = sqrt(scale.l/param.g);//graph->ticks[1]->time - graph->ticks[0]->time;
	scale.tinv = 1.0/scale.t;
	scale.at   = param.g;
	scale.vt   = scale.at*scale.t;
	scale.ft   = param.m*param.g;
	scale.pt   = scale.vt*scale.t;
	scale.pt2  = scale.pt*scale.pt;
	scale.pr   = scale.pt/scale.l;
	scale.vr   = scale.vt/scale.l;
	scale.ar   = scale.at/scale.l;
	scale.fr   = scale.ft*scale.l;
	scale.L    = scale.fr*scale.t;

	//scale.l    = 1.0;
	//scale.t    = 1.0;
	//scale.tinv = 1.0;
	//scale.at   = 1.0;
	//scale.vt   = 1.0;
	//scale.ft   = 1.0;
	//scale.pt   = 1.0;
	//scale.pr   = 1.0;
	//scale.vr   = 1.0;
	//scale.ar   = 1.0;
	//scale.fr   = 1.0;
	//scale.L    = 1.0;
}

void Centroid::CalcComAcceleration(CentroidData& d){
	vec3_t fsum = zero3;

	int nend = (int)d.ends.size();
	for(int i = 0; i < nend; i++){
		CentroidData::End& dend = d.ends[i];

		fsum += dend.force_t;
	}

	d.acc_t = (1.0/param.m)*fsum - vec3_t(0.0, 0.0, param.g);
}

void Centroid::CalcBaseAcceleration(CentroidData& d){
	vec3_t msum = zero3;

	int nend = (int)d.ends.size();
	for(int i = 0; i < nend; i++){
		CentroidData::End& dend = d.ends[i];

		msum += dend.force_r + (dend.pos_t - d.pos_t).cross(dend.force_t);
	}

	//d.acc_r = param.Iinv*msum;
	// not correct (time-varying inertia)
	d.acc_r = d.Iinv[0]*msum;
}

void Centroid::CalcWrench(CentroidData& d){
	vec3_t pc = d.pos_t;

	int nend = (int)d.ends.size();
	for(int i = 0; i < nend; i++){
		CentroidData::End& dend = d.ends[i];

		if(param.endWrenchParametrization == EndWrenchParametrization::Stiffness){
			real_t le  = dend.stiff;
			real_t le2 = le*le;
			vec2_t re  = dend.cmp;
			vec3_t me  = dend.moment;
			vec3_t pe  = dend.pos_t;

			dend.force_t = param.m*le2*(pc - pe - vec3_t(re.x(), re.y(), 0.0));
			dend.force_r = param.m*le2*me;
		}
		else{
			dend.force_t = dend.force;
			dend.force_r = dend.moment;
		}
	}
}

void Centroid::CalcStiffness(CentroidData& d){
    int nend  = (int)ends .size();
    
	if(param.endWrenchParametrization == EndWrenchParametrization::Stiffness){
		/*
		min    sum li^4
		sub.to sum li^2 (pc - pi) = ac + g

		li^2 = (pc - pi)^T mu
		(sum (pc - pi)(pc - pi)^T) mu = g
		*/
		vector<vec3_t>  pl(nend);
		mat3_t A = eps*mat3_t::Identity();
		for(int i = 0; i < nend; i++){
			pl[i] = d.pos_t - d.ends[i].pos_t;
			A += vvtrmat(pl[i], pl[i]);
		}
		mat3_t Ainv = A.inverse();

		for(int i = 0; i < nend; i++){
			real_t li = sqrt(std::max(0.0, pl[i].dot(Ainv*vec3_t(0.0, 0.0, param.g))));
			d.ends[i].stiff = li;
			d.ends[i].cmp = ends[i].cmpOffset;					
			d.ends[i].moment = zero3;
		}
	}
	else{
		for(int i = 0; i < nend; i++){
			d.ends[i].force  = vec3_t(0.0, 0.0, param.m*param.g/(real_t)nend);
			d.ends[i].moment = zero3;
		}
	}
}

int Centroid::NumSteps(){
	// count keys from phases
	int N = 0;
	for(int i = 0; i < (int)phases.size()-1; i++){
		N += phases[i].ndiv;
	}
	return N;
}

void Centroid::Reset(bool reset_first, bool reset_middle, bool reset_last){
	int nend   = (int)ends  .size();
	int N = (int)world->ticks.size()-1;

	for (int k = 0; k <= N; k++) {
		if(k == 0 && !reset_first)
			continue;
		if(0 < k && k < N && !reset_middle)
			continue;
		if(k == N && !reset_last)
			continue;

		CentroidKey* key = (CentroidKey*)traj.GetKeypoint(k);
		CentroidData& d = key->data_des;

		key->var_pos_t->val = d.pos_t;
		key->var_pos_r->val = d.pos_r;
		key->var_vel_t->val = d.vel_t;
		key->var_L    ->val = d.L;
		key->var_time ->val = d.time;

		for(int i = 0; i < nend; i++){
			CentroidKey ::End&  end  = key->ends[i];
			CentroidData::End&  dend = d.ends [i];
			
			end.var_pos_t->val = dend.pos_t;
			end.var_pos_r->val = dend.pos_r;
		}

		if(key->next){
			key->var_duration->val = d.duration;
		
			for(int i = 0; i < nend; i++){
				CentroidKey ::End&  end  = key->ends[i];
				CentroidData::End&  dend = d.ends [i];
			
				end.var_vel_t->val = dend.vel_t;
				end.var_vel_r->val = dend.vel_r;

				if(param.endWrenchParametrization == EndWrenchParametrization::Stiffness){
					end.var_stiff ->val = dend.stiff ;
					end.var_cmp[0]->val = dend.cmp[0];
					end.var_cmp[1]->val = dend.cmp[1];
					end.var_moment->val = dend.moment;
				}
				else{
					end.var_force ->val = dend.force ;
					end.var_moment->val = dend.moment;
				}
			}
		}
	}

	trajReady = false;
}

void Centroid::Shift(){
	int nend   = (int)ends  .size();
	int N = (int)world->ticks.size()-1;

	for (int k = 0; k < N; k++) {
		CentroidKey* key0 = (CentroidKey*)traj.GetKeypoint(k+0);
		CentroidKey* key1 = (CentroidKey*)traj.GetKeypoint(k+1);
		
		key0->var_pos_t->val = key1->var_pos_t->val;
		key0->var_pos_r->val = key1->var_pos_r->val;
		key0->var_vel_t->val = key1->var_vel_t->val;
		key0->var_L    ->val = key1->var_L    ->val;
		key0->var_time ->val = key1->var_time ->val;

		if(key1->next){
			key0->var_duration->val = key1->var_duration->val;
		}
		
		for(int i = 0; i < nend; i++){
			key0->ends[i].var_pos_t->val = key1->ends[i].var_pos_t->val;
			key0->ends[i].var_pos_r->val = key1->ends[i].var_pos_r->val;
			key0->ends[i].var_vel_t->val = key1->ends[i].var_vel_t->val;
			key0->ends[i].var_vel_r->val = key1->ends[i].var_vel_r->val;

			if(key1->next){
				if(param.endWrenchParametrization == EndWrenchParametrization::Stiffness){
					key0->ends[i].var_stiff ->val = key1->ends[i].var_stiff ->val;
					key0->ends[i].var_cmp[0]->val = key1->ends[i].var_cmp[0]->val;
					key0->ends[i].var_cmp[1]->val = key1->ends[i].var_cmp[1]->val;
					key0->ends[i].var_moment->val = key1->ends[i].var_moment->val;
				}
				else{
					key0->ends[i].var_force ->val = key1->ends[i].var_force ->val;
					key0->ends[i].var_moment->val = key1->ends[i].var_moment->val;
				}
			}
		}
	}

	trajReady = false;
}

void Centroid::Setup(){
	int nend  = (int)ends.size();
	int nface = (int)faces.size();
	int N = (int)world->ticks.size()-1;

	if(!phases.empty()){
		int iphase = 0;
		int idiv = 0;
		for(int k = 0; k <= N; k++){
			CentroidKey* key = (CentroidKey*)traj.GetKeypoint(k);
			CentroidData& d = key->data_des;

			key->iphase = iphase;
			key->idiv   = idiv;			

			for(int i = 0; i < nend; i++){
				d.ends[i].iface = phases[iphase].iface[i];
			}

			if(++idiv == phases[iphase].ndiv){
				idiv = 0;
				iphase++;
			}
		}	
    }

	// setup from waypoints
	if(!waypoints.empty()){
		curve_real_t          curve_time;
		curve_vec3_t          curve_pos_t;
		curve_vec3_t          curve_pos_r;
		curve_real_t          curve_weight_time;
		curve_real_t          curve_weight_duration;
		curve_vec3_t          curve_weight_pos_t;
		curve_vec3_t          curve_weight_pos_r;
		curve_vec3_t          curve_weight_vel_t;
		curve_vec3_t          curve_weight_L;

        curve_time           .type = Interpolate::LinearDiff;
		curve_pos_t          .type = Interpolate::Cubic     ;
		curve_pos_r          .type = Interpolate::Cubic     ;
		curve_weight_time    .type = Interpolate::LinearDiff;
		curve_weight_duration.type = Interpolate::LinearDiff;
		curve_weight_pos_t   .type = Interpolate::LinearDiff;
		curve_weight_pos_r   .type = Interpolate::LinearDiff;
		curve_weight_vel_t   .type = Interpolate::LinearDiff;
		curve_weight_L       .type = Interpolate::LinearDiff;

		vector<curve_vec3_t>  curve_end_pos_t;
		vector<curve_vec3_t>  curve_end_pos_r;
		vector<curve_vec3_t>  curve_end_weight_pos_t ;
		vector<curve_vec3_t>  curve_end_weight_pos_r ;
		vector<curve_vec3_t>  curve_end_weight_vel_t ;
		vector<curve_vec3_t>  curve_end_weight_vel_r ;
		vector<curve_real_t>  curve_end_weight_stiff ;
		vector<curve_vec3_t>  curve_end_weight_moment;
		vector<curve_vec2_t>  curve_end_weight_cmp   ;
	
		curve_end_pos_t        .resize(nend);
		curve_end_pos_r        .resize(nend);
		curve_end_weight_pos_t .resize(nend);
		curve_end_weight_pos_r .resize(nend);
		curve_end_weight_vel_t .resize(nend);
		curve_end_weight_vel_r .resize(nend);
		curve_end_weight_stiff .resize(nend);
		curve_end_weight_cmp   .resize(nend);
		curve_end_weight_moment.resize(nend);

        for(int j = 0; j < nend; j++){
			curve_end_pos_t[j].type = Interpolate::Cubic;
			curve_end_pos_r[j].type = Interpolate::Cubic;
			curve_end_weight_pos_t [j].type = Interpolate::LinearDiff;
			curve_end_weight_pos_r [j].type = Interpolate::LinearDiff;
			curve_end_weight_vel_t [j].type = Interpolate::LinearDiff;
			curve_end_weight_vel_r [j].type = Interpolate::LinearDiff;
			curve_end_weight_stiff [j].type = Interpolate::LinearDiff;
			curve_end_weight_cmp   [j].type = Interpolate::LinearDiff;
			curve_end_weight_moment[j].type = Interpolate::LinearDiff;
		}

		// interpolation of time and weights
		for(int k = 0; k < (int)waypoints.size(); k++){
			Waypoint& wp = waypoints[k];
			if(is_valid(wp.value.time)){
				curve_time.points.push_back(curve_real_t::point_t(k, wp.value.time, 0.0));
			}
			if(is_valid(wp.weight.time)){
				curve_weight_time.points.push_back(curve_real_t::point_t(k, wp.weight.time, 0.0));
			}
			if(is_valid(wp.weight.duration)){
				curve_weight_duration.points.push_back(curve_real_t::point_t(k, wp.weight.duration, 0.0));
			}
			if(is_valid(wp.weight.pos_t)){
				curve_weight_pos_t.points.push_back(curve_vec3_t::point_t(k, wp.weight.pos_t, vec3_t::Zero()));
			}
			if(is_valid(wp.weight.pos_r)){
				curve_weight_pos_r.points.push_back(curve_vec3_t::point_t(k, wp.weight.pos_r, vec3_t::Zero()));
			}
			if(is_valid(wp.weight.vel_t)){
				curve_weight_vel_t.points.push_back(curve_vec3_t::point_t(k, wp.weight.vel_t, vec3_t::Zero()));
			}
			if(is_valid(wp.weight.L)){
				curve_weight_L.points.push_back(curve_vec3_t::point_t(k, wp.weight.L, vec3_t::Zero()));
			}
			
			for(int j = 0; j < nend; j++){
				if(wp.ends.size() <= j)
					continue;

				if(is_valid(wp.ends[j].weight.pos_t)){
					curve_end_weight_pos_t[j].points.push_back(curve_vec3_t::point_t(k, wp.ends[j].weight.pos_t, vec3_t::Zero()));
				}
				if(is_valid(wp.ends[j].weight.pos_r)){
					curve_end_weight_pos_r[j].points.push_back(curve_vec3_t::point_t(k, wp.ends[j].weight.pos_r, vec3_t::Zero()));
				}
				if(is_valid(wp.ends[j].weight.vel_t)){
					curve_end_weight_vel_t[j].points.push_back(curve_vec3_t::point_t(k, wp.ends[j].weight.vel_t, vec3_t::Zero()));
				}
				if(is_valid(wp.ends[j].weight.vel_r)){
					curve_end_weight_vel_r[j].points.push_back(curve_vec3_t::point_t(k, wp.ends[j].weight.vel_r, vec3_t::Zero()));
				}
				if(is_valid(wp.ends[j].weight.stiff)){
					curve_end_weight_stiff[j].points.push_back(curve_real_t::point_t(k, wp.ends[j].weight.stiff, 0.0));
				}
				if(is_valid(wp.ends[j].weight.cmp)){
					curve_end_weight_cmp[j].points.push_back(curve_vec2_t::point_t(k, wp.ends[j].weight.cmp, vec2_t::Zero()));
				}
				if(is_valid(wp.ends[j].weight.moment)){
					curve_end_weight_moment[j].points.push_back(curve_vec3_t::point_t(k, wp.ends[j].weight.moment, vec3_t::Zero()));
				}
			}
		}

		// interpolation of trajectory
		for(int k = 0; k < (int)waypoints.size(); k++){
			Waypoint& wp = waypoints[k];
			real_t t = curve_time.CalcPos(k);

			if(is_valid(wp.value.pos_t)){
				curve_pos_t.points.push_back(curve_vec3_t::point_t(t, wp.value.pos_t, is_valid(wp.value.vel_t) ? wp.value.vel_t : zero3));
			}
			if(is_valid(wp.value.pos_r)){
				curve_pos_r.points.push_back(curve_vec3_t::point_t(t, wp.value.pos_r, is_valid(wp.value.vel_r) ? wp.value.vel_r : zero3));
			}

			for(int j = 0; j < nend; j++){
				if(wp.ends.size() <= j)
					continue;

				if(is_valid(wp.ends[j].value.pos_t)){
					curve_end_pos_t[j].points.push_back(curve_vec3_t::point_t(t, wp.ends[j].value.pos_t, is_valid(wp.ends[j].value.vel_t) ? wp.ends[j].value.vel_t : zero3));
				}
				if(is_valid(wp.ends[j].value.pos_r)){
					curve_end_pos_r[j].points.push_back(curve_vec3_t::point_t(t, wp.ends[j].value.pos_r, is_valid(wp.ends[j].value.vel_r) ? wp.ends[j].value.vel_r : zero3));
				}
			}
		}

		for (int k = 0; k <= N; k++) {
			CentroidKey* key = (CentroidKey*)traj.GetKeypoint(k);
			CentroidData& d = key->data_des;

			real_t s0 = (real_t)key->iphase + (real_t)(key->idiv+0)/(real_t)phases[key->iphase].ndiv;
			real_t s1 = (real_t)key->iphase + (real_t)(key->idiv+1)/(real_t)phases[key->iphase].ndiv;

			// initial setting of time
			real_t t  = curve_time.CalcPos(s0);
			real_t dt = curve_time.CalcPos(s1) - t;

			d.time     = t;
			d.duration = dt;
			d.time_weight     = curve_weight_time.CalcPos(s0);
			d.duration_weight = curve_weight_duration.CalcPos(s0);
		
			d.pos_t = curve_pos_t.CalcPos(t);
			d.vel_t = curve_pos_t.CalcVel(t);

			vec3_t rpy  = curve_pos_r.CalcPos(t);
			vec3_t rpyd = curve_pos_r.CalcVel(t);
			d.pos_r = FromRollPitchYaw(rpy);
			d.vel_r = VelocityFromRollPitchYaw(rpy, rpyd);
			d.L = param.I*d.vel_r;
		
			d.pos_t_weight = curve_weight_pos_t.CalcPos(s0);
			d.vel_t_weight = curve_weight_vel_t.CalcPos(s0);
			d.pos_r_weight = curve_weight_pos_r.CalcPos(s0);
			d.L_weight     = curve_weight_L    .CalcPos(s0);

			for(int i = 0; i < nend; i++){
				CentroidData::End& dend = d.ends[i];

				dend.pos_t = curve_end_pos_t[i].CalcPos(t);
				dend.vel_t = curve_end_pos_t[i].CalcVel(t);
				dend.pos_r = FromRollPitchYaw(curve_end_pos_r[i].CalcPos(t));
				dend.vel_r = curve_end_pos_r[i].CalcVel(t);

				dend.pos_t_weight = curve_end_weight_pos_t[i].CalcPos(s0);
				dend.vel_t_weight = curve_end_weight_vel_t[i].CalcPos(s0);
				dend.pos_r_weight = curve_end_weight_pos_r[i].CalcPos(s0);
				dend.vel_r_weight = curve_end_weight_vel_r[i].CalcPos(s0);
		
				dend.stiff_weight  = curve_end_weight_stiff [i].CalcPos(s0);
				dend.cmp_weight    = curve_end_weight_cmp   [i].CalcPos(s0);
				dend.moment_weight = curve_end_weight_moment[i].CalcPos(s0);
			}

			CalcStiffness(d);
		}
	}

    // setup using callback
	if(callback){
		for (int k = 0; k <= N; k++) {
			CentroidKey* key = (CentroidKey*)traj.GetKeypoint(k);

			real_t t = world->ticks[k]->time;

			// need to get contact state as desired state
			callback->GetDesiredState(k, t, key->data_des);
			CalcStiffness(key->data_des);
		
			if(k == 0){
				callback->GetInitialState(key->data_des);
			}

			// copy inertia info
			key->data.I    = key->data_des.I;
			key->data.Iinv = key->data_des.Iinv;
		}
	}

	for (int k = 0; k <= N; k++) {
		CentroidKey* key = (CentroidKey*)traj.GetKeypoint(k);
		CentroidData& d = key->data_des;

		if(k == 0){
			key->var_pos_t->val = d.pos_t;
			key->var_pos_r->val = d.pos_r;
			key->var_vel_t->val = d.vel_t;
			key->var_L    ->val = d.L;
			key->var_time ->val = d.time;
		}
	
		key->con_des_pos_t->desired = d.pos_t;
		key->con_des_pos_r->desired = d.pos_r;
		key->con_des_vel_t->desired = d.vel_t;
		key->con_des_L    ->desired = d.L;
		key->con_des_pos_t->weight  = d.pos_t_weight;
		key->con_des_pos_r->weight  = d.pos_r_weight;
		key->con_des_vel_t->weight  = d.vel_t_weight;
		key->con_des_L    ->weight  = d.L_weight;

		key->con_des_time    ->desired   = d.time;
		key->con_des_time    ->weight[0] = d.time_weight;

		if(key->next){
			int ndiv = (!phases.empty() ? phases[key->iphase].ndiv : 1);
			key->con_duration_range[0]->bound =  (param.durationMin/ndiv);
			key->con_duration_range[1]->bound = -(param.durationMax/ndiv);
			key->con_duration_range[0]->weight[0] = 0.1;
			key->con_duration_range[1]->weight[0] = 0.1;
			key->con_duration_range[0]->barrier_margin = 0.00001;
			key->con_duration_range[1]->barrier_margin = 0.00001;

			key->con_des_duration->desired   = d.duration;
			key->con_des_duration->weight[0] = d.duration_weight;
		}
		
		for(int i = 0; i < nend; i++){
			CentroidKey ::End&  end  = key->ends[i];
			CentroidData::End&  dend = d.ends [i];

			if(k == 0){
				end.var_pos_t->val = dend.pos_t;
				end.var_vel_t->val = (dend.iface == -1 ? dend.vel_t : zero3);
				end.var_pos_r->val = dend.pos_r;
				end.var_vel_r->val = (dend.iface == -1 ? dend.vel_r : zero3);			
			}

			for(int j = 0; j < 3; j++){
				key->ends[i].con_pos_range[j][0]->weight[0] = 1.0;
				key->ends[i].con_pos_range[j][1]->weight[0] = 1.0;
				key->ends[i].con_pos_range[j][0]->barrier_margin = 0.001;
				key->ends[i].con_pos_range[j][1]->barrier_margin = 0.001;
			}

			end.con_des_pos_t->desired = dend.pos_t;
			end.con_des_vel_t->desired = (dend.iface == -1 ? dend.vel_t : zero3);
			end.con_des_pos_r->desired = dend.pos_r;
			end.con_des_vel_r->desired = (dend.iface == -1 ? dend.vel_r : zero3);
			end.con_des_pos_t->weight  = dend.pos_t_weight;
			end.con_des_vel_t->weight  = (dend.iface == -1 ? dend.vel_t_weight : param.complWeight*one3);
			end.con_des_pos_r->weight  = dend.pos_r_weight;
			end.con_des_vel_r->weight  = (dend.iface == -1 ? dend.vel_r_weight : param.complWeight*one3);
					
			if(key->next){
				if(param.endWrenchParametrization == EndWrenchParametrization::Stiffness){
					key->ends[i].con_stiff_range->_min = 0.0;
					key->ends[i].con_stiff_range->_max = ends[i].stiffnessMax;
					key->ends[i].con_stiff_range->weight[0] = 0.1;
					key->ends[i].con_stiff_range->barrier_margin = 0.00001;

					end.con_des_stiff ->desired   = (dend.iface != -1 ? dend.stiff  : 0.0);
					end.con_des_cmp[0]->desired   = (dend.iface != -1 ? dend.cmp[0] : 0.0);
					end.con_des_cmp[1]->desired   = (dend.iface != -1 ? dend.cmp[1] : 0.0);
					end.con_des_moment->desired   = (dend.iface != -1 ? dend.moment : zero3);
					end.con_des_stiff ->weight[0] = (dend.iface != -1 ? dend.stiff_weight  : param.complWeight);
					end.con_des_cmp[0]->weight[0] = (dend.iface != -1 ? dend.cmp_weight[0] : param.complWeight);
					end.con_des_cmp[1]->weight[0] = (dend.iface != -1 ? dend.cmp_weight[1] : param.complWeight);
					end.con_des_moment->weight    = (dend.iface != -1 ? dend.moment_weight : param.complWeight*one3);
				}
				else{
					end.con_des_force ->weight = (dend.iface != -1 ? dend.force_weight  : param.complWeight*one3);
					end.con_des_moment->weight = (dend.iface != -1 ? dend.moment_weight : param.complWeight*one3);
				}
				key->ends[i].con_friction[0][0]->weight[0] = 0.1;
				key->ends[i].con_friction[0][1]->weight[0] = 0.1;
				key->ends[i].con_friction[1][0]->weight[0] = 0.1;
				key->ends[i].con_friction[1][1]->weight[0] = 0.1;
				key->ends[i].con_friction[0][0]->barrier_margin = 0.00001;
				key->ends[i].con_friction[0][1]->barrier_margin = 0.00001;
				key->ends[i].con_friction[1][0]->barrier_margin = 0.00001;
				key->ends[i].con_friction[1][1]->barrier_margin = 0.00001;
				key->ends[i].con_moment[0][0]->weight[0] = 1.0;
				key->ends[i].con_moment[0][1]->weight[0] = 1.0;
				key->ends[i].con_moment[1][0]->weight[0] = 1.0;
				key->ends[i].con_moment[1][1]->weight[0] = 1.0;
				key->ends[i].con_moment[2][0]->weight[0] = 1.0;
				key->ends[i].con_moment[2][1]->weight[0] = 1.0;
				key->ends[i].con_moment[0][0]->barrier_margin = 0.0001;
				key->ends[i].con_moment[0][1]->barrier_margin = 0.0001;
				key->ends[i].con_moment[1][0]->barrier_margin = 0.0001;
				key->ends[i].con_moment[1][1]->barrier_margin = 0.0001;
				key->ends[i].con_moment[2][0]->barrier_margin = 0.0001;
				key->ends[i].con_moment[2][1]->barrier_margin = 0.0001;
			}

			for(int j = 0; j < nface; j++){
				end.con_contact[j]->weight[0] = (dend.iface == j ? param.complWeight : 0.0);
			}
		}
	}
}

void Centroid::Init() {
	Model::Init();

	param.Iinv = param.I.inverse();

    int nend  = (int)ends .size();
    int nface = (int)faces.size();
	int N     = (int)world->ticks.size()-1;

	for (int k = 0; k <= N; k++) {
		CentroidKey* key = (CentroidKey*)traj.GetKeypoint(k);

		key->var_pos_r  ->locked = !param.enableRotation;
		key->var_L      ->locked = !param.enableRotation;

		for(int i = 0; i < nend; i++){
			key->ends[i].var_pos_r->locked = !param.enableRotation || ends[i].lockOri;
			key->ends[i].var_vel_r->locked = !param.enableRotation || ends[i].lockOri;
		}

		if(key->next){
			for(int i = 0; i < nend; i++){
				if(param.endWrenchParametrization == EndWrenchParametrization::Stiffness){
					for(int j = 0; j < 2; j++){
						key->ends[i].var_cmp[j]->locked = ends[i].lockCmp;
					}
				}
				key->ends[i].var_moment->locked = ends[i].lockMoment;
            }
		}
	}

	// call prepare here so that initial trajectory is visualized properly
    Prepare();

	trajReady = false;
}

void Centroid::Prepare() {
	trajReady = false;
	traj.Update();

	#pragma omp parallel for if(world->solver->param.parallelize)
	for(int k = 0; k < traj.size(); k++){
		traj[k]->Prepare();
	}

	#pragma omp parallel for if(world->solver->param.parallelize)
	for(int k = 0; k < traj.size(); k++){
		((CentroidKey*)traj.GetKeypoint(k))->Prepare2();
	}
}

void Centroid::PrepareStep() {
	#pragma omp parallel for if(world->solver->param.parallelize)
	for(int k = 0; k < traj.size(); k++){
		traj[k]->PrepareStep();
	}
}

void Centroid::Finish(){
	Model::Finish();
}

void Centroid::CalcState(real_t t, CentroidData& d){
	if(traj.empty())
		return;

	KeyPair      kp = traj.GetSegment(t);
	CentroidKey* k0 = (CentroidKey*)kp.first;
	CentroidKey* k1 = (CentroidKey*)kp.second;

	CalcState(t, k0->data, k1->data, d);
}

void Centroid::CalcState(real_t t, const vector<CentroidData>& d_array, CentroidData& d){
	int N = (int)d_array.size()-1;

	if(t < d_array[0].time){
		CalcState(t, d_array[0], d_array[0], d);
		return;
	}

	for(int k = 0; k < N; k++){
		if(t < d_array[k+1].time){
			CalcState(t, d_array[k+0], d_array[k+1], d);
			return;
		}
	}

	CalcState(t, d_array[N], d_array[N], d);
}

void Centroid::CalcState(real_t t, const CentroidData& d0, const CentroidData& d1, CentroidData& d){
	int nend = (int)ends.size();

	d.ends.resize(nend);

	d.time = t;
	d.duration = d0.duration;

	if(d1.time == d0.time){
        d.pos_t = d0.pos_t;
        d.vel_t = d0.vel_t;
		d.acc_t = zero3;

		d.pos_r = d0.pos_r;
		d.L     = d0.L;
        if(d.Llocal.empty())
			d.Llocal.push_back(d0.Llocal[0]);

		for(int i = 0; i < nend; i++){
			d.ends[i].pos_t = d0.ends[i].pos_t;
			d.ends[i].pos_r = d0.ends[i].pos_r;
			d.ends[i].vel_t = d0.ends[i].vel_t;
			d.ends[i].vel_r = d0.ends[i].vel_r;
		}
	}
	else{
		real_t t0  = d0.time;
		real_t t1  = d1.time;
		real_t dt  = t - t0;
		real_t dt2 = dt*dt;
		real_t dt3 = dt*dt2;

		if(param.endWrenchParametrization == EndWrenchParametrization::Stiffness){
			real_t Ct = cosh(d0.lbar*(dt));
			real_t St = sinh(d0.lbar*(dt));
	
			d.pos_t = d0.pbar + d0.rbar + Ct*(d0.pos_t - (d0.pbar + d0.rbar)) + (St/d0.lbar)*d0.vel_t;
			d.vel_t = d0.lbar * St*(d0.pos_t - (d0.pbar + d0.rbar)) + Ct*d0.vel_t;
			d.acc_t = (d0.lbar*d0.lbar)*(d.pos_t - (d0.pbar + d0.rbar));
		
			d.L = d0.L + param.m*( (d.vel_t - d0.vel_t).cross(d0.rbar) + (t - t0)*d0.etabar);
		}
		else{
			d.acc_t = d0.fsum/param.m - vec3_t(0.0, 0.0, param.g);
			d.vel_t = d0.vel_t + d.acc_t*dt;
			d.pos_t = d0.pos_t + d0.vel_t*dt + d.acc_t*(dt*dt/2.0);

			d.L = d0.L + d0.etasum*dt + d0.fsum.cross(d0.pos_t*dt + d0.vel_t*(dt2/2.0) - (dt3/6.0)*vec3_t(0.0, 0.0, param.g));
		}

		int ndiv    = param.rotationResolution;
		real_t dtau = d0.duration/(real_t)ndiv;
		int idiv = std::min(std::max(0, (int)((t - t0)/dtau)), ndiv-1);
		
		d.pos_r = interpolate_slerp_diff(
			t,
			t0 + (idiv+0)*dtau, d0.q_rhs[idiv+0],
			t0 + (idiv+1)*dtau, d0.q_rhs[idiv+1]);
		
		d.vel_r = d0.Iinv[idiv]*(d.L - d0.Llocal[idiv]);

        if(d.Llocal.empty())
			d.Llocal.push_back(d0.Llocal[idiv]);

		for(int i = 0; i < nend; i++){
			if(d0.ends[i].state == ContactState::Free){
				if(param.endInterpolation == EndInterpolation::Polynomial){
					real_t t0, t1, _t;
					t0 = (d0.ends[i].t_lift + param.swingLiftMargin < d0.time ? d0.time : d0.ends[i].t_lift + param.swingLiftMargin);
					t1 = (d1.ends[i].t_land - param.swingLandMargin > d1.time ? d1.time : d1.ends[i].t_land - param.swingLandMargin);
					_t = std::min(std::max(t0, t), t1);

					d.ends[i].pos_t = interpolate_pos_cubic(_t, t0, d0.ends[i].pos_t, d0.ends[i].vel_t_ave, t1, d1.ends[i].pos_t, d1.ends[i].vel_t_ave);
					d.ends[i].vel_t = interpolate_pos_cubic(_t, t0, d0.ends[i].pos_t, d0.ends[i].vel_t_ave, t1, d1.ends[i].pos_t, d1.ends[i].vel_t_ave);
					
					t0 = (d0.ends[i].t_lift + param.swingTurnMargin < d0.time ? d0.time : d0.ends[i].t_lift + param.swingTurnMargin);
					t1 = (d1.ends[i].t_land - param.swingTurnMargin > d1.time ? d1.time : d1.ends[i].t_land - param.swingTurnMargin);
					_t = std::min(std::max(t0, t), t1);

					vec3_t angle = interpolate_pos_cubic(_t, 
						t0, ToRollPitchYaw(d0.ends[i].pos_r), d0.ends[i].vel_r_ave, 
						t1, ToRollPitchYaw(d1.ends[i].pos_r), d1.ends[i].vel_r_ave);
					vec3_t angled = interpolate_vel_cubic(_t, 
						t0, ToRollPitchYaw(d0.ends[i].pos_r), d0.ends[i].vel_r_ave, 
						t1, ToRollPitchYaw(d1.ends[i].pos_r), d1.ends[i].vel_r_ave);
					d.ends[i].pos_r = FromRollPitchYaw(angle);
					d.ends[i].vel_r = VelocityFromRollPitchYaw(angle, angled);
					//d.ends[i].pos_r = InterpolateOri   (t, d0.time, d0.ends[i].pos_r, vec3_t(), d1.time, d1.ends[i].pos_r, vec3_t(), Interpolate::SlerpDiff);
					//d.ends[i].vel_r = InterpolateAngvel(t, d0.time, d0.ends[i].pos_r, vec3_t(), d1.time, d1.ends[i].pos_r, vec3_t(), Interpolate::SlerpDiff);
				}
				else{
				    const real_t _2pi = 2.0*_pi;
				    const real_t tmargin = 0.05;
				    real_t tau = (d1.ends[i].t_land - d0.ends[i].t_lift) - 2.0*tmargin;
				    real_t s   = std::min(std::max(0.0, (t - (d0.ends[i].t_lift + tmargin))/tau), 1.0);
				    real_t ch  = (s - sin(_2pi*s)/_2pi);
				    real_t chd = ((1.0 - cos(_2pi*s))/tau);
				    real_t cv  = (1 - cos(_2pi*s))/2.0;
				    real_t cvd = (_2pi*sin(_2pi*s)/(2.0*tau));

				    real_t sw = param.swingHeight;
				    vec3_t nx = d1.ends[i].pos_t_land - d0.ends[i].pos_t_lift;
				    vec3_t ny(0.0, 1.0, 0.0);
				    vec3_t nz(0.0, 0.0, 1.0);
				    real_t nxnorm = nx.norm();

				    if(nxnorm > eps){
					    nx = nx/nxnorm;
					    ny = nz.cross(nx);
					    nz = nx.cross(ny);
				    }
			
				    if(param.endInterpolation == EndInterpolation::CycloidLocal){
					    vec3_t pe = d0.ends[i].pos_t_lift_local + ch *(d1.ends[i].pos_t_land_local - d0.ends[i].pos_t_lift_local) + (cv *sw)*nz;
					    vec3_t ve =                               chd*(d1.ends[i].pos_t_land_local - d0.ends[i].pos_t_lift_local) + (cvd*sw)*nz;

					    Eigen::AngleAxisd qrel(d0.ends[i].pos_r_lift_local.conjugate()*d1.ends[i].pos_r_land_local);
					    vec3_t axis   = qrel.axis ();
					    real_t theta  = qrel.angle();
					    if(theta > _pi)
						    theta -= 2*_pi;
				
					    quat_t qe = d0.ends[i].pos_r_lift*rot_quat((ch*theta)*axis);
					    vec3_t we = d0.ends[i].pos_r_lift*((chd*theta)*axis);

					    d.ends[i].pos_t = d.pos_t + d.pos_r*pe;
					    d.ends[i].pos_r = d.pos_r*qe;
					    d.ends[i].vel_t = d.vel_t + d.pos_r*ve + d.vel_r.cross(d.pos_r*pe);
					    d.ends[i].vel_r = d.vel_r + d.pos_r*we;
				    }
				    else{
					    /*
						(p1 - p0)*(p0 + s(p1 - p0) - pc) = 0
						s = (pc - p0)*(p1 - p0)/||p1 - p0||^2
						*/
						// calc side sway based on distance between straight swing line and support foot
						vec2_t p0(d0.ends[i].pos_t_lift.x(), d0.ends[i].pos_t_lift.y());
						vec2_t p1(d1.ends[i].pos_t_land.x(), d1.ends[i].pos_t_land.y());
						vec2_t psup(d0.ends[!i].pos_t.x(), d0.ends[!i].pos_t.y());
						vec2_t swh = zero2;
						/*real_t s = (psup - p0)*(p1 - p0)/(p1 - p0).square();
						if(0.0 < s && s < 1.0){
							vec2_t psup_proj = p0 + s*(p1 - p0);
							const real_t th = 0.1;
							real_t r = (psup_proj - psup).norm();
							if(r < th){
								swh = ((th - r)/r)*(psup_proj - psup);
							}
						}*/
						
						vec3_t pe = d0.ends[i].pos_t_lift + ch *(d1.ends[i].pos_t_land - d0.ends[i].pos_t_lift) + (cv *sw)*nz + (cv )*vec3_t(swh.x(), swh.y(), 0.0);
						vec3_t ve =                         chd*(d1.ends[i].pos_t_land - d0.ends[i].pos_t_lift) + (cvd*sw)*nz + (cvd)*vec3_t(swh.x(), swh.y(), 0.0);

						Eigen::AngleAxisd qrel(d0.ends[i].pos_r_lift.conjugate()*d1.ends[i].pos_r_land);
					    vec3_t axis   = qrel.axis ();
					    real_t theta  = qrel.angle();
					    if(theta > _pi)
						    theta -= 2*_pi;
								
						quat_t qe = d0.ends[i].pos_r_lift*rot_quat((ch*theta)*axis);
						vec3_t we = d0.ends[i].pos_r_lift*((chd*theta)*axis);

						d.ends[i].pos_t = pe;
						d.ends[i].pos_r = qe;
						d.ends[i].vel_t = ve;
						d.ends[i].vel_r = we;
				    }
                }
			}
			// surface contact
			else if(d0.ends[i].state == ContactState::Surface && d1.ends[i].state == ContactState::Surface){
				d.ends[i].pos_t = d0.ends[i].pos_t;
				d.ends[i].pos_r = d0.ends[i].pos_r;
				d.ends[i].vel_t = d0.ends[i].vel_t;
				d.ends[i].vel_r = d0.ends[i].vel_r;
			}
			// other contact
			else{
				d.ends[i].pos_t = interpolate_pos_linear_diff(t, d0.time, d0.ends[i].pos_t, d1.time, d1.ends[i].pos_t);
				d.ends[i].vel_t = interpolate_vel_linear_diff(   d0.time, d0.ends[i].pos_t, d1.time, d1.ends[i].pos_t);
				d.ends[i].pos_r = interpolate_slerp_diff     (t, d0.time, d0.ends[i].pos_r, d1.time, d1.ends[i].pos_r);
				d.ends[i].vel_r = interpolate_angvel_diff    (   d0.time, d0.ends[i].pos_r, d1.time, d1.ends[i].pos_r);
			}
		}
	}

	for(int i = 0; i < nend; i++){
		d.ends[i].iface  = d0.ends[i].iface;
		if(param.endWrenchParametrization == EndWrenchParametrization::Stiffness){
			d.ends[i].stiff  = d0.ends[i].stiff;
			d.ends[i].cmp    = d0.ends[i].cmp;
			d.ends[i].moment = d0.ends[i].moment;
		}
		else{
			d.ends[i].force  = d0.ends[i].force;
			d.ends[i].moment = d0.ends[i].moment;
		}
	}

	CalcWrench(d);
}	

void Centroid::CalcTrajectory() {
	real_t tf = traj.back()->tick->time;
	real_t dt = 0.01;

	trajectory.clear();
	for (real_t t = 0.0; t <= tf; t += dt) {
		CentroidData d;
		d.Init(this);
		CalcState(t, d);
		trajectory.push_back(d);
	}

	trajReady = true;
}

void Centroid::Draw(render::Canvas* canvas, render::Config* conf) {
	Model::Draw(canvas, conf);

	if (!trajReady)
		CalcTrajectory();

	if (trajectory.empty())
		return;

	// pos
	if (conf->Set(canvas, render::Item::CentroidPos, this)) {
		canvas->BeginLayer("centroid_pos", true);
		canvas->BeginPath();
		canvas->MoveTo(trajectory[0].pos_t);
		for(int i = 1; i < trajectory.size(); i++){
			canvas->LineTo(trajectory[i].pos_t);
		}
		canvas->EndPath();
		canvas->EndLayer();
	}

	stringstream ss;

	// end
	if (conf->Set(canvas, render::Item::CentroidEndTraj, this)) {
		for(int i = 0; i < ends.size(); i++){
			ss.str("");
			ss << i;
            canvas->SetLineColor(i % 2 ? "blue" : "magenta");
			canvas->BeginLayer("centroid_end" + ss.str(), true);
			canvas->BeginPath();
			canvas->MoveTo(trajectory[0].ends[i].pos_t);
			for (int k = 1; k < trajectory.size(); k++) {
				canvas->LineTo(trajectory[k].ends[i].pos_t);
			}
			canvas->EndPath();
			canvas->EndLayer();
		}
	}
    /*
    if(conf->Set(canvas, Render::Item::CentroidFace, this)){
	    for(int i = 0; i < faces.size(); i++){
			Face& f = faces[i];
			ss.str("");
			ss << i;
			canvas->BeginLayer("centroid_face" + ss.str(), true);
		    canvas->BeginPath();
		    canvas->MoveTo(f.hull->vertices[0]);
		    canvas->LineTo(f.hull->vertices[1]);
		    canvas->LineTo(f.hull->vertices[2]);
		    canvas->LineTo(f.hull->vertices[3]);
		    canvas->LineTo(f.hull->vertices[0]);
		    canvas->EndPath();
		    canvas->EndLayer();
	    }
    }
    */
}

void Centroid::CreateSnapshot(real_t t){
	snapshot.Init(this);
	CalcState(t, snapshot);
}

void Centroid::DrawSnapshot(render::Canvas* canvas, render::Config* conf) {
	if (conf->Set(canvas, render::Item::CentroidEnd, this)) {
		canvas->BeginPath();
		canvas->MoveTo(snapshot.pos_t);
		canvas->LineTo(snapshot.pos_t + snapshot.pos_r*vec3_t(0.0, 0.0, 0.3));
		canvas->EndPath();

		for(int i = 0; i < ends.size(); i++){
			canvas->BeginLayer("centroid_end_snapshot", true);
			canvas->SetLineColor("black");
            canvas->SetLineWidth(1.0f);

			// line connecting com and end
		    canvas->BeginPath();
			canvas->MoveTo(snapshot.pos_t);
			canvas->LineTo(snapshot.pos_t + snapshot.pos_r*ends[i].basePos);
			canvas->LineTo(snapshot.ends[i].pos_t);
			canvas->EndPath();

            // line indicating force
            canvas->SetLineColor("red");
            canvas->SetLineWidth(1.0f);
			canvas->BeginPath();
			canvas->MoveTo(snapshot.ends[i].pos_t);
			canvas->LineTo(snapshot.ends[i].pos_t + 0.001*snapshot.ends[i].force_t);
			canvas->EndPath();
			canvas->EndLayer();
		}

        // end rectangle
        for(int i = 0; i < ends.size(); i++){
            vec3_t vtx[4];
            vtx[0] = vec3_t(ends[i].copMin.x(), ends[i].copMin.y(), 0.0);
            vtx[1] = vec3_t(ends[i].copMin.x(), ends[i].copMax.y(), 0.0);
            vtx[2] = vec3_t(ends[i].copMax.x(), ends[i].copMax.y(), 0.0);
            vtx[3] = vec3_t(ends[i].copMax.x(), ends[i].copMin.y(), 0.0);
            canvas->SetLineColor("green");
            canvas->SetLineWidth(/*snapshot.ends[i].contact ? 2.0f : */1.0f);
			canvas->BeginPath();
			canvas->MoveTo(snapshot.ends[i].pos_t + snapshot.ends[i].pos_r*vtx[0]);
			canvas->LineTo(snapshot.ends[i].pos_t + snapshot.ends[i].pos_r*vtx[1]);
			canvas->LineTo(snapshot.ends[i].pos_t + snapshot.ends[i].pos_r*vtx[2]);
			canvas->LineTo(snapshot.ends[i].pos_t + snapshot.ends[i].pos_r*vtx[3]);
			canvas->LineTo(snapshot.ends[i].pos_t + snapshot.ends[i].pos_r*vtx[0]);
			canvas->EndPath();
		}
	}	
}

///////////////////////////////////////////////////////////////////////////////////////////////////

CentroidCon::CentroidCon(Solver* solver, int _dim, int _tag, string _name, CentroidKey* _obj, real_t _scale):
	Constraint(solver, _dim, ID(_tag, _obj->model, _name), Constraint::Type::Equality, _scale) {
	obj[0] = _obj;
	obj[1] = (CentroidKey*)_obj->next;
}

CentroidComCon::CentroidComCon(Solver* solver, int _dim, int _tag, string _name, CentroidKey* _obj, real_t _scale):
	CentroidCon(solver, _dim, _tag, _name, _obj, _scale) {

}

CentroidPosConT::CentroidPosConT(Solver* solver, string _name, CentroidKey* _obj, real_t _scale):
	CentroidComCon(solver, 3, ConTag::CentroidPosT, _name, _obj, _scale) {

	AddSLink (obj[1]->var_pos_t);
	AddSLink (obj[0]->var_pos_t);
	AddSLink (obj[0]->var_vel_t);
	AddC3Link(obj[0]->var_duration);
	
	int nend = (int)obj[0]->ends.size();
	for(int i = 0; i < nend; i++){
		if(obj[0]->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			AddSLink (obj[0]->ends[i].var_pos_t);
			AddSLink (obj[0]->ends[i].var_vel_t);
			AddC3Link(obj[0]->ends[i].var_stiff);
			AddC3Link(obj[0]->ends[i].var_cmp[0]);
			AddC3Link(obj[0]->ends[i].var_cmp[1]);
		}
		else{
			AddSLink (obj[0]->ends[i].var_force);
		}
	}
}

CentroidVelConT::CentroidVelConT(Solver* solver, string _name, CentroidKey* _obj, real_t _scale):
	CentroidComCon(solver, 3, ConTag::CentroidVelT, _name, _obj, _scale) {

	AddSLink (obj[1]->var_vel_t);
	AddSLink (obj[0]->var_pos_t);
	AddSLink (obj[0]->var_vel_t);
	AddC3Link(obj[0]->var_duration);

	int nend = (int)obj[0]->ends.size();
	for(int i = 0; i < nend; i++){
		if(obj[0]->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			AddSLink (obj[0]->ends[i].var_pos_t);
			AddSLink (obj[0]->ends[i].var_vel_t);
			AddC3Link(obj[0]->ends[i].var_stiff);
			AddC3Link(obj[0]->ends[i].var_cmp[0]);
			AddC3Link(obj[0]->ends[i].var_cmp[1]);
		}
		else{
			AddSLink (obj[0]->ends[i].var_force);
		}
	}
}

CentroidPosConR::CentroidPosConR(Solver* solver, string _name, CentroidKey* _obj, real_t _scale):
	CentroidCon(solver, 3, ConTag::CentroidPosR, _name, _obj, _scale) {

	AddSLink (obj[1]->var_pos_r);
	AddM3Link(obj[0]->var_pos_r);
	AddM3Link(obj[0]->var_L    );
	AddM3Link(obj[0]->var_pos_t);
	AddM3Link(obj[0]->var_vel_t);
	AddC3Link(obj[0]->var_duration);

	int nend = (int)obj[0]->ends.size();
	for(int i = 0; i < nend; i++){
		AddM3Link(obj[0]->ends[i].var_pos_t );
		AddM3Link(obj[0]->ends[i].var_vel_t );
		if(obj[0]->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			AddC3Link(obj[0]->ends[i].var_stiff );
			AddC3Link(obj[0]->ends[i].var_cmp[0]);
			AddC3Link(obj[0]->ends[i].var_cmp[1]);
			AddM3Link(obj[0]->ends[i].var_moment);
		}
		else{
			AddM3Link(obj[0]->ends[i].var_force );
			AddM3Link(obj[0]->ends[i].var_moment);
		}
	}
}

CentroidLCon::CentroidLCon(Solver* solver, string _name, CentroidKey* _obj, real_t _scale):
	CentroidCon(solver, 3, ConTag::CentroidMomentum, _name, _obj, _scale) {

	AddSLink (obj[1]->var_L);
	AddSLink (obj[0]->var_L);
	AddM3Link(obj[0]->var_pos_t);
	AddM3Link(obj[0]->var_vel_t);
	AddC3Link(obj[0]->var_duration);	

	int nend = (int)obj[0]->ends.size();
	for(int i = 0; i < nend; i++){
		AddM3Link(obj[0]->ends[i].var_pos_t );
		AddM3Link(obj[0]->ends[i].var_vel_t);
		if(obj[0]->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			AddC3Link(obj[0]->ends[i].var_stiff );
			AddC3Link(obj[0]->ends[i].var_cmp[0]);
			AddC3Link(obj[0]->ends[i].var_cmp[1]);
			AddSLink (obj[0]->ends[i].var_moment);
		}
		else{
			AddM3Link(obj[0]->ends[i].var_force );
			AddSLink (obj[0]->ends[i].var_moment);
		}
	}
}

CentroidTimeCon::CentroidTimeCon(Solver* solver, string _name, CentroidKey* _obj, real_t _scale):
	CentroidCon(solver, 1, ConTag::CentroidTime, _name, _obj, _scale) {
	
	AddSLink(obj[1]->var_time);
	AddSLink(obj[0]->var_time);
	AddSLink(obj[0]->var_duration);
}

CentroidEndPosConT::CentroidEndPosConT(Solver* solver, string _name, CentroidKey* _obj, int _iend, real_t _scale):
	CentroidCon(solver, 3, ConTag::CentroidEndPosT, _name, _obj, _scale) {
	iend  = _iend;

	AddSLink (obj[1]->ends[iend].var_pos_t);
	AddSLink (obj[0]->ends[iend].var_pos_t);
	AddSLink (obj[0]->ends[iend].var_vel_t);
	AddC3Link(obj[0]->var_duration);
}

CentroidEndPosConR::CentroidEndPosConR(Solver* solver, string _name, CentroidKey* _obj, int _iend, real_t _scale):
	CentroidCon(solver, 3, ConTag::CentroidEndPosR, _name, _obj, _scale) {
	iend  = _iend;

	AddSLink (obj[1]->ends[iend].var_pos_r);
	AddM3Link(obj[0]->ends[iend].var_pos_r);
	AddM3Link(obj[0]->ends[iend].var_vel_r);
	AddC3Link(obj[0]->var_duration);
}

CentroidDurationRangeCon::CentroidDurationRangeCon(Solver* solver, string _name, CentroidKey* _obj, real_t _dir, real_t _scale):
	Constraint(solver, 1, ID(ConTag::CentroidDurationRange, _obj->model, _name), Constraint::Type::InequalityBarrier, _scale){
	obj  = _obj;
	dir  = _dir;

	AddSLink(obj->var_duration);
}

CentroidEndPosRangeCon::CentroidEndPosRangeCon(Solver* solver, string _name, CentroidKey* _obj, int _iend, int _idx, real_t _dir, real_t _scale):
	Constraint(solver, 1, ID(ConTag::CentroidEndPosRange, _obj->model, _name), Constraint::Type::InequalityBarrier, _scale){
	obj  = _obj;
	iend = _iend;
	idx  = _idx;
	dir  = _dir;

	AddR3Link(obj->var_pos_t);
	AddR3Link(obj->var_pos_r);
	AddR3Link(obj->ends[iend].var_pos_t);
}

CentroidEndContactCon::CentroidEndContactCon(Solver* solver, string _name, CentroidKey* _obj, int _iend, int _iface, real_t _scale):
	Constraint(solver, 1, ID(ConTag::CentroidEndContact, _obj->model, _name), Constraint::Type::Equality, _scale){
	obj   = _obj;
	iend  = _iend;
    iface = _iface;
    face  = &obj->cen->faces[_iface];

	AddR3Link(obj->ends[iend].var_pos_t);
}

CentroidEndFrictionCon::CentroidEndFrictionCon(Solver* solver, string _name, CentroidKey* _obj, int _iend, int _dir, int _side, real_t _scale):
	Constraint(solver, 1, ID(ConTag::CentroidEndFriction, _obj->model, _name), Constraint::Type::InequalityBarrier, _scale){
	obj  = _obj;
	iend = _iend;
    dir  = _dir;
	side = _side;
    
	if(obj->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
		AddR3Link(obj->var_pos_t);
		AddR3Link(obj->ends[iend].var_pos_t);
		AddR3Link(obj->ends[iend].var_vel_t);
		AddSLink (obj->ends[iend].var_cmp[0]);
		AddSLink (obj->ends[iend].var_cmp[1]);
	}
	else{
		AddR3Link(obj->ends[iend].var_force);
	}
}

CentroidEndMomentCon::CentroidEndMomentCon(Solver* solver, string _name, CentroidKey* _obj, int _iend, int _dir, int _side, real_t _scale):
	Constraint(solver, 1, ID(ConTag::CentroidEndMomentRange, _obj->model, _name), Constraint::Type::InequalityBarrier, _scale){
	obj   = _obj;
	iend  = _iend;
    dir   = _dir;
	side  = _side;
    
	if(obj->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
		AddR3Link(obj->var_pos_t);
		AddR3Link(obj->ends[iend].var_pos_t );
		AddSLink (obj->ends[iend].var_cmp[0]);
		AddSLink (obj->ends[iend].var_cmp[1]);
		AddR3Link(obj->ends[iend].var_moment);
	}
	else{
		AddR3Link(obj->ends[iend].var_force );
		AddR3Link(obj->ends[iend].var_moment);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CentroidCon::Prepare(){
}

void CentroidComCon::Prepare(){
	CentroidCon::Prepare();

}

void CentroidPosConT::Prepare(){
	CentroidComCon::Prepare();
}

void CentroidVelConT::Prepare(){
	CentroidComCon::Prepare();
}

void CentroidPosConR::Prepare(){
	CentroidCon::Prepare();
}

void CentroidLCon::Prepare(){
	CentroidCon::Prepare();
}

void CentroidTimeCon::Prepare(){
	CentroidCon::Prepare();

	t     = obj[0]->data.time;
	t_lhs = obj[1]->data.time;
	t_rhs = t + obj[0]->data.duration;
}

void CentroidEndPosConT::Prepare(){
	CentroidCon::Prepare();

	pe     = obj[0]->data.ends[iend].pos_t;
	pe_lhs = obj[1]->data.ends[iend].pos_t;
	ve     = obj[0]->data.ends[iend].vel_t;
	pe_rhs = pe + obj[0]->data.duration*ve;
}

void CentroidEndPosConR::Prepare(){
	CentroidCon::Prepare();

	qe     = obj[0]->data.ends[iend].pos_r;
	qe_lhs = obj[1]->data.ends[iend].pos_r;

	we    = obj[0]->data.ends[iend].vel_r;
	omega = obj[0]->data.duration*we;
	q_omega = rot_quat(omega);
	R_omega = q_omega.toRotationMatrix();
	A_omega = rot_jacobian(omega);
	qe_rhs  = q_omega*qe;
}

void CentroidEndPosRangeCon::Prepare(){
	p       = obj->data.pos_t;
	q       = obj->data.pos_r;
	pe      = obj->data.ends[iend].pos_t;
	pbase   = obj->cen->ends[iend].basePos;
	
	eta = zero3;
	eta[idx] = dir;
	eta_abs = q*eta;

	if(dir == 1.0)
		 bound =  obj->cen->ends[iend].posMin[idx];
	else bound = -obj->cen->ends[iend].posMax[idx];

}

void CentroidEndContactCon::Prepare(){
    pe = obj->data.ends[iend].pos_t;
    pf = vec3_t(pe.x(), pe.y(), 0.0);
    nf = face->normal;
}

void CentroidEndFrictionCon::Prepare(){
	/*
	 f = p - pi - ri   // ignore  m*li^2
	 mu*fz - sqrt(fx^2 + fy^2) >= 0
	*/
	if(obj->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness)
		 f = obj->data.pos_t - (obj->pi[iend] + obj->ri[iend]);
	else f = obj->fi[iend];

	mu = obj->cen->param.mu;
	qi = obj->data.ends[iend].pos_r;
	nx = qi*ex;
	ny = qi*ey;
	nz = qi*ez;
	ft = (dir == 0 ? nx : ny).dot(f);
	fz = nz.dot(f);
	df = (side == 0 ? -1.0 : 1.0)*(dir == 0 ? nx : ny) + mu*nz;
		
	/*
	ft = vec2_t(f.x, f.y);
	ftnorm = ft.norm();
	if(ftnorm < eps)
		 ftn = vec2_t();
	else ftn = ft/ftnorm;

	mu = obj->cen->param.mu;

	d = vec3_t(-ftn.x, -ftn.y, mu);
	*/
}

void CentroidEndMomentCon::Prepare(){
	if(obj->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
		f   = obj->data.pos_t - (obj->pi[iend] + obj->ri[iend]);
		eta = obj->etai[iend];
	}
	else{
		f   = obj->fi[iend];
		eta = obj->etai[iend];
	}
	
	qi   = obj->data.ends[iend].pos_r;
	nx   = qi*ex;
	ny   = qi*ey;
	nz   = qi*ez;
	etax = nx.dot(eta);
	etay = ny.dot(eta);
	etaz = nz.dot(eta);
	fz   = std::max(nz.dot(f), 0.0);
	cmin = obj->cen->ends[iend].copMin;
	cmax = obj->cen->ends[iend].copMax;

	if(dir == 0){
		df   = (side == 0 ? -cmin.x() : cmax.x())*nz;
		deta = (side == 0 ? -1.0 :  1.0)*ny;
	}
	if(dir == 1){
		df   = (side == 0 ? -cmin.y() : cmax.y())*nz;
		deta = (side == 0 ?  1.0 : -1.0)*nx;
	}
	if(dir == 2){
		df   = (side == 0 ? -cmin.z() : cmax.z())*nz;
		deta = (side == 0 ?  1.0 : -1.0)*nz;
	}

	/*Rc  = mat3_t();

	if(dir.x > 0.0)
		 bound = vec3_t(obj->cen->ends[iend].copMin.y, -obj->cen->ends[iend].copMax.x, -obj->cen->param.mu);
	else bound = vec3_t(obj->cen->ends[iend].copMax.y, -obj->cen->ends[iend].copMin.x,  obj->cen->param.mu);
	*/
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CentroidPosConT::CalcCoef(){
	Prepare();

	int i = 0;
	dynamic_cast<SLink *>(links[i++])->SetCoef( 1.0    );
	dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->p_p  );
	dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->p_v  );
	dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->p_tau);
	
	int nend = (int)obj[0]->ends.size();
	for(int iend = 0; iend < nend; iend++){
		if(obj[0]->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->p_pi[iend]);
			dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->p_vi[iend]);
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->p_li[iend]);
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->p_ri[iend]*ex);
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->p_ri[iend]*ey);
		}
		else{
			dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->p_fi[iend]);
		}
	}
}

void CentroidVelConT::CalcCoef(){
	Prepare();

	int i = 0;
	dynamic_cast<SLink *>(links[i++])->SetCoef( 1.0    );
	dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->v_p  .back());
	dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->v_v  .back());
	dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->v_tau.back());

	int nend = (int)obj[0]->ends.size();
	int ndiv = obj[0]->cen->param.rotationResolution;
	for(int iend = 0; iend < nend; iend++){
		if(obj[0]->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->v_pi[iend][ndiv]);
			dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->v_vi[iend][ndiv]);
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->v_li[iend][ndiv]);
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->v_ri[iend][ndiv]*ex);
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->v_ri[iend][ndiv]*ey);
		}
		else{
			dynamic_cast<SLink *>(links[i++])->SetCoef(-obj[0]->v_fi[iend][ndiv]);
		}
	}
}

void CentroidPosConR::CalcCoef(){
	Prepare();

	int i = 0;
	dynamic_cast<SLink* >(links[i++])->SetCoef( 1.0);
	dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->q_q  );
	dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->q_L  );
	dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->q_p  );
	dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->q_v  );
	dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->q_tau);

	int nend = (int)obj[0]->ends.size();
	for(int iend = 0; iend < nend; iend++){
		dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->q_pi  [iend]);
		dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->q_vi  [iend]);
		if(obj[0]->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->q_li  [iend]);
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->q_ri  [iend].col(0));
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->q_ri  [iend].col(1));
			dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->q_etai[iend]);
		}
		else{
			dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->q_fi  [iend]);
			dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->q_etai[iend]);
		}
	}	
}

void CentroidLCon::CalcCoef(){
	Prepare();

	int nend = (int)obj[0]->ends.size();
	int ndiv = obj[0]->cen->param.rotationResolution;

	int i = 0;
	dynamic_cast<SLink* >(links[i++])->SetCoef( 1.0);
	dynamic_cast<SLink* >(links[i++])->SetCoef(-obj[0]->L_L  );
	dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->L_p  [ndiv]);
	dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->L_v  [ndiv]);
	dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->L_tau[ndiv]);
	
	for(int iend = 0; iend < nend; iend++){
		dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->L_pi  [iend][ndiv]);
		dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->L_vi  [iend][ndiv]);
		if(obj[0]->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->L_li  [iend][ndiv]);
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->L_ri  [iend][ndiv].col(0));
			dynamic_cast<C3Link*>(links[i++])->SetCoef(-obj[0]->L_ri  [iend][ndiv].col(1));
			dynamic_cast<SLink* >(links[i++])->SetCoef(-obj[0]->L_etai[iend][ndiv]);
		}
		else{
			dynamic_cast<M3Link*>(links[i++])->SetCoef(-obj[0]->L_fi  [iend][ndiv]);
			dynamic_cast<SLink* >(links[i++])->SetCoef(-obj[0]->L_etai[iend][ndiv]);
		}
	}	
}

void CentroidTimeCon::CalcCoef(){
	Prepare();

	((SLink*)links[0])->SetCoef( 1.0);
	((SLink*)links[1])->SetCoef(-1.0);
	((SLink*)links[2])->SetCoef(-1.0);
}

void CentroidEndPosConT::CalcCoef(){
	Prepare();

	int i = 0;
	((SLink *)links[i++])->SetCoef( 1.0);
	((SLink *)links[i++])->SetCoef(-1.0);
	((SLink *)links[i++])->SetCoef(-obj[0]->data.duration);
	((C3Link*)links[i++])->SetCoef(-ve);
}

void CentroidEndPosConR::CalcCoef(){
	Prepare();

	int i = 0;
	((SLink *)links[i++])->SetCoef( 1.0);
	((M3Link*)links[i++])->SetCoef(-R_omega);
	((M3Link*)links[i++])->SetCoef(-A_omega*obj[0]->data.duration);
	((C3Link*)links[i++])->SetCoef(-A_omega*we);
}

void CentroidDurationRangeCon::CalcCoef(){
	((SLink*)links[0])->SetCoef(dir);
}

void CentroidEndPosRangeCon::CalcCoef(){
	Prepare();

	/* 
	(q*dir)^T *(pe - (p + q*pbase)) >= bound
	(q*dir)^T *(pe - p) - dir^T*pbase >= bound
	
	(pe - p)^T*(q*dir)
	  (pe - p)^T*(Omega % dir_abs)
	= (pe - p)^T*(dir_abs^xT)*Omega
	= (dir_abs % (pe - p))^T *Omega
	*/

	((R3Link*)links[0])->SetCoef(-eta_abs);
	((R3Link*)links[1])->SetCoef( eta_abs.cross(pe - p));
	((R3Link*)links[2])->SetCoef( eta_abs);
}

void CentroidEndContactCon::CalcCoef(){
	Prepare();

	((R3Link*)links[0])->SetCoef(nf);
}

void CentroidEndFrictionCon::CalcCoef(){
	Prepare();

	int idx = 0;
	if(obj->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
		dynamic_cast<R3Link*>(links[idx++])->SetCoef( df);
		dynamic_cast<R3Link*>(links[idx++])->SetCoef(-df);
		dynamic_cast<R3Link*>(links[idx++])->SetCoef(-df*(obj->data.duration/2.0));
		dynamic_cast<SLink* >(links[idx++])->SetCoef(-df.x());
		dynamic_cast<SLink* >(links[idx++])->SetCoef(-df.y());
	}
	else{
		dynamic_cast<R3Link*>(links[idx++])->SetCoef( df);
	}
}

void CentroidEndMomentCon::CalcCoef(){
	Prepare();
	
	// f   = l2*(p - pi - ri)
	// eta = l2*etahat

	// -my - (cop_min.x)*fz >= 0
	//  my + (cop_max.x)*fz >= 0
	//  mx - (cop_min.y)*fz >= 0
	// -mx + (cop_max.y)*fz >= 0
	//  mz - (cop_min.z)*fz >= 0
	// -mz + (cop_max.z)*fz >= 0
	
	int idx = 0;
	if(obj->cen->param.endWrenchParametrization == Centroid::EndWrenchParametrization::Stiffness){
		dynamic_cast<R3Link*>(links[idx++])->SetCoef( df);
		dynamic_cast<R3Link*>(links[idx++])->SetCoef(-df);
		dynamic_cast<SLink *>(links[idx++])->SetCoef(-df.x());
		dynamic_cast<SLink *>(links[idx++])->SetCoef(-df.y());
		dynamic_cast<R3Link*>(links[idx++])->SetCoef( deta);

		//dynamic_cast<R3Link*>(links[idx++])->SetCoef(-(dir*bound)*Rc.col(2));
		//dynamic_cast<R3Link*>(links[idx++])->SetCoef( (dir*bound)*Rc.col(2));
		//dynamic_cast<R3Link*>(links[idx++])->SetCoef( (Rc*dir) % eta - (dir*bound)*(Rc.col(2) % f) );
		//dynamic_cast<SLink *>(links[idx++])->SetCoef( (dir*bound)*Rc.col(2)[0]);
		//dynamic_cast<SLink *>(links[idx++])->SetCoef( (dir*bound)*Rc.col(2)[1]);
		//dynamic_cast<R3Link*>(links[idx++])->SetCoef( Rc*dir);
	}
	else{
		dynamic_cast<R3Link*>(links[idx++])->SetCoef(df  );
		dynamic_cast<R3Link*>(links[idx++])->SetCoef(deta);
		
		//dynamic_cast<R3Link*>(links[idx++])->SetCoef(vec3_t(0.0, 0.0, -(dir*bound)*f.z));
		//dynamic_cast<R3Link*>(links[idx++])->SetCoef(dir);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void CentroidPosConT::CalcDeviation(){
	y = obj[1]->data.pos_t - obj[0]->data.p_rhs;
}

void CentroidPosConR::CalcDeviation(){
	Eigen::AngleAxisd qerror(obj[0]->data.q_rhs.back().conjugate()*obj[1]->data.pos_r);
	vec3_t axis   = qerror.axis ();
	real_t theta  = qerror.angle();
	if(theta > _pi)
		theta -= 2*_pi;

	y = (1.0/2.0)*(obj[0]->data.q_rhs.back()*(theta*axis) + obj[1]->data.pos_r*(theta*axis));
}

void CentroidVelConT::CalcDeviation(){
	y = obj[1]->data.vel_t - obj[0]->data.v_rhs.back();
}
void CentroidLCon::CalcDeviation(){
	y = obj[1]->data.L - obj[0]->data.L_rhs.back();
}
void CentroidTimeCon::CalcDeviation(){
	y[0] = t_lhs - t_rhs;
}
void CentroidEndPosConT::CalcDeviation(){
	y = pe_lhs - pe_rhs;
}

void CentroidEndPosConR::CalcDeviation(){
	Eigen::AngleAxisd qerror(qe_rhs.conjugate()*qe_lhs);
	vec3_t axis   = qerror.axis ();
	real_t theta  = qerror.angle();
	if(theta > _pi)
		theta -= 2*_pi;

	y = (1.0/2.0)*(qe_rhs*(theta*axis) + qe_lhs*(theta*axis));
}

void CentroidDurationRangeCon::CalcDeviation(){
    y[0] = dir*obj->data.duration - bound;

	// set activity if penalty mode
	active = (y[0] < 0.0);
}

void CentroidEndPosRangeCon::CalcDeviation(){
    y[0] = eta_abs.dot(pe - (p + q*pbase)) - bound;
	
	// set activity if penalty mode
	if(y[0] < 0.0){
		active = true;
	}
	else{
		active = false;
	}
}

void CentroidEndContactCon::CalcDeviation(){
	y[0] = nf.dot(pe - pf) - obj->cen->param.contactMargin;
}

void CentroidEndFrictionCon::CalcDeviation(){
	y[0] = df.dot(f);
	//y[0] = mu*f.z - ftnorm;
	active = (obj->data_des.ends[iend].iface != -1 && y[0] < 0.0);
}

void CentroidEndMomentCon::CalcDeviation(){
	y[0] = df.dot(f) + deta.dot(eta);
	//y[0] = dir*(eta - bound*f.z);
}

}
