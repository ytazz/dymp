#include <solver.h>
#include <variable.h>
#include <constraint.h>
#include <util.h>

#define USE_MKL

#if defined USE_MKL
# ifdef _WIN32
#  include <mkl_lapacke.h>
# else
//#  include <lapacke.h>
#  include <mkl_lapacke.h>
# endif
#endif

namespace dymp{;

///////////////////////////////////////////////////////////////////////////////////////////////////

Solver::SubInput* Solver::Input::Find(Variable* var){
	for(auto&& subin : subinput){
		if(subin->var == var)
			return subin.get();
	}
	return 0;
}

Solver::SubState* Solver::State::Find(Variable* var){
	for(auto&& subst : substate){
		if(subst->var == var)
			return subst.get();
	}
	return 0;
}

Solver::SubState* Solver::AddStateVar(Variable* var, int k){
	if(state.size() <= k)
		state.resize(k+1);

	if(!state[k])
		state[k] = unique_ptr<State>(new State());

	SubState* subst = new SubState();
	subst->var = var;
	state[k]->substate.push_back(unique_ptr<SubState>(subst));

    return subst;
}

Solver::SubInput* Solver::AddInputVar(Variable* var, int k){
	if(input.size() <= k)
		input.resize(k+1);

	if(!input[k])
		input[k] = unique_ptr<Input>(new Input());

	SubInput* subin = new SubInput();
	subin->var = var;
	input[k]->subinput.push_back(unique_ptr<SubInput>(subin));

    return subin;
}

Solver::SubTransition* Solver::AddTransitionCon(Constraint* con, int k){
	if(transition.size() <= k)
		transition.resize(k+1);

	if(!transition[k])
		transition[k] = unique_ptr<Transition>(new Transition());

	auto&& st_x0 = state[k+0];
	auto&& st_x1 = state[k+1];
	auto&& in_u  = input[k+0];

	SubTransition* subtr = new SubTransition();
	subtr->con = con;
	for(Link* l : con->links){
		SubState* subst;
		SubInput* subin;

		// variable that belongs to st_x1 must be unique
		subst = st_x1->Find(l->var);
		if(subst){
			subtr->x1 = subst;
		}
		subst = st_x0->Find(l->var);
		if(subst){
			SubStateLink xl;
			xl.x    = subst;
			xl.link = l;
			subtr->x0.push_back(xl);
		}
		subin = in_u->Find(l->var);
		if(subin){
			SubInputLink ul;
			ul.u    = subin;
			ul.link = l;
			subtr->u.push_back(ul);
		}
	}

	transition[k]->subtran.push_back(unique_ptr<SubTransition>(subtr));

    return subtr;
}

Solver::SubCost* Solver::AddCostCon (Constraint* con, int k){
	if(cost.size() <= k)
		cost.resize(k+1);

	if(!cost[k])
		cost[k] = unique_ptr<Cost>(new Cost());

	State* st_x = (k < state.size() ? (State*)state[k].get() : (State*)0);
	Input* in_u = (k < input.size() ? (Input*)input[k].get() : (Input*)0);

	SubCost* subcost = new SubCost();
	subcost->con = con;

	for(Link* l : con->links){
		SubState* subst;
		SubInput* subin;

		subst = (st_x ? st_x->Find(l->var) : (SubState*)0);
		if(subst){
			SubStateLink xl;
			xl.x    = subst;
			xl.link = l;
			subcost->x.push_back(xl);
		}
		subin = (in_u ? in_u->Find(l->var) : (SubInput*)0);
		if(subin){
			SubInputLink ul;
			ul.u    = subin;
			ul.link = l;
			subcost->u.push_back(ul);
		}
	}

	cost[k]->subcost.push_back(unique_ptr<SubCost>(subcost));

    return subcost;
}

void Solver::InitDDP(){
	N = (int)state.size()-1;

	for(int k = 0; k < state.size(); k++){
        auto&& st = state[k];
		st->dim = 0;

		for(auto&& subst : st->substate){
			if(subst->var->locked)
				continue;

			subst->index = st->dim;
			st->dim += subst->var->nelem;
		}
	}

	for(int k = 0; k < input.size(); k++){
        auto&& in = input[k];
		in->dim = 0;

		for(auto&& subin : in->subinput){
			if(subin->var->locked)
				continue;

			subin->index = in->dim;
			in->dim += subin->var->nelem;
		}
	}

	for(int k = 0; k <= N; k++){
		if(!cost[k])
			continue;

		cost[k]->dim = 0;

		for(auto&& subcost : cost[k]->subcost){     
			subcost->index = -1;

			if(!subcost->con->enabled)
				continue;

			subcost->xbegin = state[k]->dim;
			subcost->xend   = 0;
			subcost->ubegin = 0;
			subcost->uend   = 0;
			for(SubStateLink& x : subcost->x){
				if(x.x->var->locked)
					continue;
				subcost->xbegin = std::min(subcost->xbegin, x.x->index);
				subcost->xend   = std::max(subcost->xend  , x.x->index + x.x->var->nelem);
			}
			if(k < N){
				subcost->ubegin = input[k]->dim;
				subcost->uend   = 0;
				for(SubInputLink& u : subcost->u){
					if(u.u->var->locked)
						continue;
					subcost->ubegin = std::min(subcost->ubegin, u.u->index);
					subcost->uend   = std::max(subcost->uend  , u.u->index + u.u->var->nelem);
				}
			}    
			if(subcost->xbegin >= subcost->xend && subcost->ubegin >= subcost->uend)
				continue;

			subcost->index = cost[k]->dim;
			cost[k]->dim += subcost->con->nelem;			        
		}

		cost[k]->y.Allocate(cost[k]->dim);
		cost[k]->b.Allocate(cost[k]->dim);
		cost[k]->Ax.Allocate(cost[k]->dim, state[k]->dim);
		if(k < N){
			cost[k]->Au.Allocate(cost[k]->dim, input[k]->dim);
		}
	}

	dx       .resize(N+1);
	du       .resize(N);
	fx       .resize(N);
	fu       .resize(N);
	fcor     .resize(N);
	L        .resize(N+1);
	Lx       .resize(N+1);
	Lxx      .resize(N+1);
	Lu       .resize(N);
	Luu      .resize(N);
	Lux      .resize(N);
	Q        .resize(N);
	Qx       .resize(N);
	Qu       .resize(N);
	Qxx      .resize(N);
	Quu      .resize(N);
	Quuinv   .resize(N);
	Quuinv_Qu.resize(N);
	Qux      .resize(N);
	V        .resize(N+1);
	Vx       .resize(N+1);
	Vxx      .resize(N+1);
	dV       .resize(N+1);
	dVx      .resize(N+1);
	dVxx     .resize(N+1);
	Vxx_fcor .resize(N);
	Vxx_fx   .resize(N);
	Vxx_fu   .resize(N);
	Vx_plus_Vxx_fcor.resize(N);
	Quuinv_Qux      .resize(N);
	Qu_plus_Qux_dx  .resize(N);

	for(int k = 0; k <= N; k++){
		int nx  = state[k]->dim;

		dx  [k].Allocate(nx);
		Lx  [k].Allocate(nx);
		Lxx [k].Allocate(nx, nx);
		Vx  [k].Allocate(nx);
		Vxx [k].Allocate(nx, nx);
		dVx [k].Allocate(nx);
		dVxx[k].Allocate(nx, nx);
		if(k == 0){
			Vxxinv.Allocate(nx, nx);
		}

		if(k < N){
			int nu  = input[k]->dim;
			int nx1 = state[k+1]->dim;

			du[k].Allocate(nu);

			fx  [k].Allocate(nx1, nx);
			fu  [k].Allocate(nx1, nu);
			fcor[k].Allocate(nx1    );
			Lu       [k].Allocate(nu      );
			Luu      [k].Allocate(nu , nu );
			Lux      [k].Allocate(nu , nx );
			Qx       [k].Allocate(nx      );
			Qu       [k].Allocate(nu      );
			Qxx      [k].Allocate(nx , nx );
			Quu      [k].Allocate(nu , nu );
			Quuinv   [k].Allocate(nu , nu );
			Quuinv_Qu[k].Allocate(nu      );
			Qux      [k].Allocate(nu , nx );
	        Vxx_fcor [k].Allocate(nx1     );
	        Vxx_fx   [k].Allocate(nx1, nx );
	        Vxx_fu   [k].Allocate(nx1, nu );
			Vx_plus_Vxx_fcor[k].Allocate(nx1);
			Quuinv_Qux      [k].Allocate(nu, nx);
			Qu_plus_Qux_dx  [k].Allocate(nu);
		}
	}

}

void Solver::ClearDDP(){
	state     .clear();
	input     .clear();
	transition.clear();
	cost      .clear();
	
	N = 0;
	dx       .clear();
	du       .clear();
	fx       .clear();
	fu       .clear();
	fcor     .clear();
	L        .clear();
	Lx       .clear();
	Lxx      .clear();
	Lu       .clear();
	Luu      .clear();
	Lux      .clear();
	Q        .clear();
	Qx       .clear();
	Qu       .clear();
	Qxx      .clear();
	Quu      .clear();
	Qux      .clear();
	Quuinv   .clear();
	Quuinv_Qu.clear();
    V        .clear();
	Vx       .clear();
	Vxx      .clear();

}

void Solver::CalcTransitionDDP(){
	timer3.CountUS();

	vec3_t one(1.0, 1.0, 1.0);

//#pragma omp parallel for if(param.parallelize)
	for(int k = 0; k < N; k++){
		auto&& tr = transition[k];
		
		mat_clear(fx  [k]);
		mat_clear(fu  [k]);
        vec_clear(fcor[k]);
		for(auto&& subtr : tr->subtran){
			if(!subtr->con->enabled)
				continue;
			if(!subtr->con->active)
				continue;
            if(subtr->x1->var->locked)
                continue;
			
			int i0 = subtr->con->index;
			int n  = subtr->con->nelem;

			// set fx
			int ix0 = 0;
			for(SubStateLink& x0 : subtr->x0){
				if(x0.x->var->locked)
					continue;

				int m = x0.x->var->nelem;
                x0.link->RegisterCoef(fx[k].SubMatrix(subtr->x1->index, x0.x->index, n, m), -one);
				
				ix0++;
			}
            
            // set fu
			for(SubInputLink& u : subtr->u){
				if(u.u->var->locked)
					continue;

				int m = u.u->var->nelem;
                u.link->RegisterCoef(fu[k].SubMatrix(subtr->x1->index, u.u->index, n, m), -one);
			}
            
            // set correction term
			subtr->con->RegisterCorrection(fcor[k].SubVector(subtr->x1->index, n), one);
		}
	}
	status.timeTrans = timer3.CountUS();

	//DSTR << " tf: " << tf << endl;
	FILE* file;
	int k = 5;
	/*
	file = fopen("fx.csv", "w");
	for(int i = 0; i < fx[k].m; i++){
		for(int j = 0; j < fx[k].n; j++){
			fprintf(file, "%f, ", fx[k](i,j));
		}
		fprintf(file, "\n");
	}
	fclose(file);

	file = fopen("fu.csv", "w");
	for(int i = 0; i < fu[k].m; i++){
		for(int j = 0; j < fu[k].n; j++){
			fprintf(file, "%f, ", fu[k](i,j));
		}
		fprintf(file, "\n");
	}
	fclose(file);
	
	file = fopen("fcor.csv", "w");
	for(int i = 0; i < fcor[k].n; i++){
		fprintf(file, "%f\n", fcor[k](i));
	}
	fclose(file);
	*/
	
}

void Solver::CalcCostDDP(){
	for(int k = 0; k <= N; k++){
		L  [k] = 0.0;
	}

	timer3.CountUS();
	// calculate state cost
//#pragma omp parallel for if(param.parallelize)
	for(int k = 0; k <= N; k++){
		if(!cost[k])
			continue;

		vec_clear(cost[k]->y);

		for(auto&& subcost : cost[k]->subcost){     
			if(subcost->index == -1)
				continue;
			if(!subcost->con->active)
				continue;
			
			int n  = subcost->con->nelem;
			subcost->con->RegisterDeviation(cost[k]->y.SubVector(subcost->index, n));

			// sum up L
			if( subcost->con->type == Constraint::Type::Equality ||
				subcost->con->type == Constraint::Type::InequalityPenalty ){
				for(int i = 0; i < n; i++){
					L[k] += 0.5 * square(subcost->con->weight[i]*cost[k]->y(subcost->index + i));
				}
			}
			if( subcost->con->type == Constraint::Type::InequalityBarrier ){
				// assume n = 1
				L[k] += -square(subcost->con->weight[0])*log(std::min(std::max(subcost->con->barrier_margin, cost[k]->y(subcost->index)), 1.0));
			}

		}
	}

	status.timeCost = timer3.CountUS();
	//DSTR << " tL: " << tL << endl;
}

void Solver::CalcCostGradientDDP(){
	for(int k = 0; k <= N; k++){
		vec_clear(Lx [k]);
		mat_clear(Lxx[k]);
	}
	for(int k = 0; k < N; k++){
		vec_clear(Lu [k]);
		mat_clear(Luu[k]);
		mat_clear(Lux[k]);
	}

	timer3.CountUS();
	// calculate state cost
//#pragma omp parallel for if(param.parallelize)
	for(int k = 0; k <= N; k++){
		if(!cost[k])
			continue;

		timer4.CountUS();
		vec_clear(cost[k]->b );
		mat_clear(cost[k]->Ax);
		mat_clear(cost[k]->Au);

		for(auto&& subcost : cost[k]->subcost){     
			if(subcost->index == -1)
				continue;
			if(!subcost->con->active)
				continue;
			
			int n  = subcost->con->nelem;

			// dynamic weight scaling
			real_t tmp;
			if( subcost->con->type == Constraint::Type::InequalityBarrier ){
				if(subcost->con->barrier_margin < cost[k]->y(subcost->index) && cost[k]->y(subcost->index) < 1.0)
					 tmp = 1.0/cost[k]->y(subcost->index);
				else tmp = 0.0;
			}
			else tmp = 1.0;

			// calc b
			if( subcost->con->type == Constraint::Type::Equality ||
				subcost->con->type == Constraint::Type::InequalityPenalty ){
				for(int i = 0; i < n; i++){
					cost[k]->b(subcost->index + i) = subcost->con->weight[i]*cost[k]->y(subcost->index + i);
				}
			}
			if( subcost->con->type == Constraint::Type::InequalityBarrier ){
				cost[k]->b(subcost->index) = -subcost->con->weight[0];
			}

			// calc Ax
			for(SubStateLink& x : subcost->x){
				if(x.x->var->locked)
					continue;

				int m = x.x->var->nelem;
				x.link->RegisterCoef(cost[k]->Ax.SubMatrix(subcost->index, x.x->index, n, m), tmp*subcost->con->weight);
			}

			// calc Au
			if(k < N){
				for(SubInputLink& u : subcost->u){
					if(u.u->var->locked)
						continue;

					int m  = u.u->var->nelem;
					u.link->RegisterCoef(cost[k]->Au.SubMatrix(subcost->index, u.u->index, n, m), tmp*subcost->con->weight);
				}
			}
		}
		int T1 = timer4.CountUS();

		timer4.CountUS();
		
		for(auto&& subcost : cost[k]->subcost){     
			if(subcost->index == -1)
				continue;
			if(!subcost->con->active)
				continue;
			
			int n  = subcost->con->nelem;
			int nx = subcost->xend - subcost->xbegin;

			if(nx > 0){
				mattr_vec_mul(
					cost[k]->Ax.SubMatrix(subcost->index, subcost->xbegin, n, nx), 
					cost[k]->b .SubVector(subcost->index, n),
					Lx[k].SubVector(subcost->xbegin, nx), 1.0, 1.0);
				mattr_mat_mul(
					cost[k]->Ax.SubMatrix(subcost->index, subcost->xbegin, n, nx), 
					cost[k]->Ax.SubMatrix(subcost->index, subcost->xbegin, n, nx), 
					Lxx[k].SubMatrix(subcost->xbegin, subcost->xbegin, nx, nx), 1.0, 1.0);
			}
			
			if(k < N){
				int nu = subcost->uend - subcost->ubegin;
				if(nu > 0){	
					mattr_vec_mul(
						cost[k]->Au.SubMatrix(subcost->index, subcost->ubegin, n, nu),
						cost[k]->b .SubVector(subcost->index, n),
						Lu[k].SubVector(subcost->ubegin, nu), 1.0, 1.0);
					mattr_mat_mul(
						cost[k]->Au.SubMatrix(subcost->index, subcost->ubegin, n, nu),
						cost[k]->Au.SubMatrix(subcost->index, subcost->ubegin, n, nu),
						Luu[k].SubMatrix(subcost->ubegin, subcost->ubegin, nu, nu), 1.0, 1.0);
					if(nx > 0){
						mattr_mat_mul(
							cost[k]->Au.SubMatrix(subcost->index, subcost->ubegin, n, nu),
							cost[k]->Ax.SubMatrix(subcost->index, subcost->xbegin, n, nx),
							Lux[k].SubMatrix(subcost->ubegin, subcost->xbegin, nu, nx), 1.0, 1.0);
					}
				}
			}
		}
		
		/*
		mattr_vec_mul(cost[k]->Ax, cost[k]->b , Lx[k], 1.0, 0.0);
		mattr_mat_mul(cost[k]->Ax, cost[k]->Ax, Lxx[k], 1.0, 0.0);
		
		if(k < N){
			mattr_vec_mul(cost[k]->Au, cost[k]->b , Lu [k], 1.0, 0.0);
			mattr_mat_mul(cost[k]->Au, cost[k]->Au, Luu[k], 1.0, 0.0);
			mattr_mat_mul(cost[k]->Au, cost[k]->Ax, Lux[k], 1.0, 0.0);
		}
		*/
		int T2 = timer4.CountUS();

		//DSTR << "costgrad: T1: " << T1 << " T2: " << T2 << endl;

        // weights
		for(auto&& subst : state[k]->substate){
            if(subst->var->locked)
                continue;

            int j0 = subst->index;
            int m  = subst->var->nelem;
            for(int j = 0; j < m; j++){
                real_t wj = subst->var->weight[j];

				Lxx[k](j0+j, j0+j) += wj*wj;
            }
        }
        if(k < N){
            for(auto&& subin : input[k]->subinput){
                if(subin->var->locked)
                    continue;

                int j0 = subin->index;
                int m  = subin->var->nelem;
                for(int j = 0; j < m; j++){
                    real_t wj = subin->var->weight[j];

					Luu[k](j0+j, j0+j) += wj*wj;
                }
            }
        }
	}
	status.timeCostGrad = timer3.CountUS();

	//DSTR << " tLgrad: " << tLgrad << endl;
	/*
	FILE* file = fopen("Luu.csv", "w");
	int k = 5;
	for(int i = 0; i < Luu[k].m; i++){
		for(int j = 0; j < Luu[k].n; j++){
			fprintf(file, "%f, ", Luu[k](i,j));
		}
		fprintf(file, "\n");
	}
	fclose(file);

	file = fopen("Lux.csv", "w");
	for(int i = 0; i < Lux[k].m; i++){
		for(int j = 0; j < Lux[k].n; j++){
			fprintf(file, "%f, ", Lux[k](i,j));
		}
		fprintf(file, "\n");
	}
	fclose(file);
	*/
}

void PrintSparsity(int k, const Matrix& m){
	FILE* file;
	char filename [256];
	sprintf(filename, "quu_%d.csv", k);
	file = fopen(filename, "w");

	for(int i = 0; i < m.m; i++){
		for(int j = 0; j < m.n; j++){
			fprintf(file, "%d", (std::abs(m(i,j)) < eps ? 0 : 1));
			if(j != m.n-1)
				fprintf(file, ", ");
		}
		fprintf(file, "\n");
	}

	fclose(file);
}

void Solver::BackwardDDP(){
	int tback1, tback2, tback3;

 	V  [N] = L  [N];
	vec_copy(Lx [N], Vx [N]);
	mat_copy(Lxx[N], Vxx[N]);
	for(int i = 0; i < Vxx[N].m; i++)
		Vxx[N](i,i) += param.stateRegularization;

	for(int k = N-1; k >= 0; k--){
		timer4.CountUS();
		mat_vec_mul(Vxx[k+1], fcor[k], Vxx_fcor[k], 1.0, 0.0);

		if(state[k]->dim != 0)
			mat_mat_mul(Vxx[k+1], fx  [k], Vxx_fx  [k], 1.0, 0.0);

		if(input[k]->dim != 0)
			mat_mat_mul(Vxx[k+1], fu  [k], Vxx_fu  [k], 1.0, 0.0);

        vec_copy(Vx[k+1]    , Vx_plus_Vxx_fcor[k]);
        vec_add (Vxx_fcor[k], Vx_plus_Vxx_fcor[k]);

        Q[k] = L[k] + V[k+1] + vec_dot(Vx[k+1], fcor[k]) + (1.0/2.0)*vec_dot(fcor[k], Vxx_fcor[k]);
        
        if(state[k]->dim != 0){
			vec_copy(Lx[k], Qx[k]);
			mattr_vec_mul(fx[k], Vx_plus_Vxx_fcor[k], Qx[k], 1.0, 1.0);
		}

        if(input[k]->dim != 0){
			vec_copy(Lu[k], Qu[k]);
			mattr_vec_mul(fu[k], Vx_plus_Vxx_fcor[k], Qu[k], 1.0, 1.0);
		}

        if(state[k]->dim != 0){
			mat_copy(Lxx[k], Qxx[k]);
	        mattr_mat_mul(fx[k], Vxx_fx[k], Qxx[k], 1.0, 1.0);
			
			for(int i = 0; i < Qxx[k].m; i++)
				Qxx[k](i,i) += param.stateRegularization;

		}

        if(input[k]->dim != 0){
			mat_copy(Luu[k], Quu[k]);
			mattr_mat_mul(fu[k], Vxx_fu[k], Quu[k], 1.0, 1.0);
		}

        if(input[k]->dim != 0 && state[k]->dim != 0){
			mat_copy(Lux[k], Qux[k]);
	        mattr_mat_mul(fu[k], Vxx_fx[k], Qux[k], 1.0, 1.0);
		}

		//Q  [k] = L  [k] + V[k+1] + Vx[k+1]*f_cor[k] + (1.0/2.0)*((Vxx[k+1]*f_cor[k])*f_cor[k]);
    	//Qx [k] = Lx [k] + fx[k].trans()*(Vx [k+1] + Vxx[k+1]*f_cor[k]);
    	//Qu [k] = Lu [k] + fu[k].trans()*(Vx [k+1] + Vxx[k+1]*f_cor[k]);
    	//Qxx[k] = Lxx[k] + fx[k].trans()*Vxx[k+1]*fx[k];
    	//Quu[k] = Luu[k] + fu[k].trans()*Vxx[k+1]*fu[k];
    	//Qux[k] = Lux[k] + fu[k].trans()*Vxx[k+1]*fx[k];
    	
		for(int i = 0; i < Quu[k].m; i++)
		    Quu[k](i,i) += param.regularization;

		//PrintSparsity(k, Quu[k]);
    	tback1 = timer4.CountUS();

		// input dimension could be zero
		if(Quu[k].m > 0){
			timer4.CountUS();
            mat_inv_pd(Quu[k], Quuinv[k]);
			tback2 = timer4.CountUS();
        
			timer4.CountUS();
			symmat_vec_mul(Quuinv[k], Qu[k], Quuinv_Qu[k], 1.0, 0.0);
			//Quuinv_Qu = Quuinv*Qu;

			V[k] = Q[k] - (1.0/2.0)*vec_dot(Qu[k], Quuinv_Qu[k]);
        
			vec_copy(Qx[k], Vx[k]);
			if(state[k]->dim != 0)
				mattr_vec_mul(Qux[k], Quuinv_Qu[k], Vx[k], -1.0, 1.0);
	    
			mat_copy(Qxx[k], Vxx[k]);
			if(state[k]->dim != 0){
				symmat_mat_mul(Quuinv[k], Qux[k], Quuinv_Qux[k], 1.0, 0.0);
				mattr_mat_mul(Qux[k], Quuinv_Qux[k], Vxx[k], -1.0, 1.0);
			}
			tback3 = timer4.CountUS();
			//mat_inv_sym(Quu[k], Quuinv[k]);
			//
			//Quuinv_Qu[k] = Quuinv[k]*Qu[k];
			//
			//V  [k] = Q  [k] - (1.0/2.0)*(Qu[k]*Quuinv_Qu[k]);
			//Vx [k] = Qx [k] - Qux[k].trans()*Quuinv_Qu [k];
			//Vxx[k] = Qxx[k] - Qux[k].trans()*Quuinv[k]*Qux[k];
		}
		else if(state[k]->dim != 0){
			V  [k] = Q  [k];
			vec_copy(Qx [k], Vx [k]);
			mat_copy(Qxx[k], Vxx[k]);
		}
				
        // enforce symmetry of Uxx
	    for(int i = 1; i < Vxx[k].m; i++) for(int j = 0; j < i; j++)
		    Vxx[k](i,j) = Vxx[k](j,i);

		//int tback = timer2.CountUS();
		//DSTR << "back " << k << " : " << " nx: " << state[k]->dim << " nu: " << input[k]->dim
		//	 << " T1: " << tback1
		//	 << " T2: " << tback2
		//	 << " T3: " << tback3
		//	 << endl;
		//DSTR << "Vk: " << k << " " << V[k] << endl;
		/*
		char filename[256];
		sprintf(filename, "Quu%d.csv", k);
		FILE* file = fopen(filename, "w");
		for(int i = 0; i < Quu[k].m; i++){
			for(int j = 0; j < Quu[k].n; j++){
				fprintf(file, "%f, ", Quu[k](i,j));
			}
			fprintf(file, "\n");
		}
		fclose(file);
		*/
	}
}

void Solver::ForwardDDP(real_t alpha){
    // if the dimension of x0 is not zero, dx0 is also optimized
	//if(state[0]->dim == 0){
	if(param.fixInitialState){
 		vec_clear(dx[0]);
	}
	else{
		mat_inv_pd(Vxx[0], Vxxinv);
		symmat_vec_mul(Vxxinv, Vx[0], dx[0], -alpha, 0.0);
    
		//dx[0] = -Vxx[0].inv()*Vx[0];
	}

    for(int k = 0; k < N; k++){
		timer4.CountUS();
		if(k == 0 && param.fixInitialInput){
			vec_clear(du[0]);
		}
		else if(input[k]->dim != 0){
			vec_copy(Qu[k], Qu_plus_Qux_dx[k]);
			for(int i = 0; i < Qu_plus_Qux_dx[k].n; i++)
				Qu_plus_Qux_dx[k](i) *= alpha;

			if(state[k]->dim != 0)
				mat_vec_mul(Qux[k], dx[k], Qu_plus_Qux_dx[k], 1.0, 1.0);

			symmat_vec_mul(Quuinv[k], Qu_plus_Qux_dx[k], du[k], -1.0, 0.0);
		}

        vec_copy(fcor[k], dx[k+1]);

		if(state[k]->dim != 0)
			mat_vec_mul(fx  [k], dx[k], dx[k+1], 1.0, 1.0);

		if(input[k]->dim != 0)
			mat_vec_mul(fu  [k], du[k], dx[k+1], 1.0, 1.0);

		//du[k] = -Quuinv[k]*(Qu[k] + Qux[k]*dx[k]);
		//dx[k+1] = fx[k]*dx[k] + fu[k]*du[k] + f_cor[k];

		int tfor = timer4.CountUS();
		//DSTR << "dx[k]: " << dx[k] << endl;
		//DSTR << "for  " << k << " : " << tfor << endl;
	}
}

void Solver::CalcDirectionDDP(){
    CalcTransitionDDP();
	CalcCostDDP();
	CalcCostGradientDDP();

	timer3.CountUS();
	BackwardDDP();
	status.timeBack = timer3.CountUS();
}

real_t Solver::CalcObjectiveDDP(){
	CalcCostDDP();

	real_t Lsum = 0.0;
	for(int k = 0; k <= N; k++)
		Lsum += L[k];

	return Lsum;
}

void Solver::ModifyVariablesDDP(real_t alpha){
    ForwardDDP(alpha);
    
	for(int k = 0; k <= N; k++){
		for(auto&& subst : state[k]->substate){
			if(subst->var->locked)
				continue;

			int j0 = subst->index;
			for(int j = 0; j < subst->var->nelem; j++){
				subst->var->dx[j] = dx[k](j0+j);
			}
		}
	}
	for(int k = 0; k < N; k++){
		for(auto&& subin : input[k]->subinput){
			if(subin->var->locked)
				continue;

			int j0 = subin->index;
			for(int j = 0; j < subin->var->nelem; j++){
				subin->var->dx[j] = du[k](j0+j);
			}
		}	
	}
}

}
