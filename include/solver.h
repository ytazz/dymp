#pragma once

#include <types.h>
#include <id.h>
#include <blas.h>
#include <timer.h>

namespace dimp3{;

/** Solver
    - solver for optimization (constraint error minimization) and optimal control problems
 **/

class Variable;
class Constraint;
class Link;

typedef std::vector< Variable* >                    Variables;
typedef std::vector< std::unique_ptr<Variable> >    VariableRefs;
typedef std::vector< Constraint* >                  Constraints;
typedef std::vector< std::unique_ptr<Constraint> >  ConstraintRefs;
typedef std::vector< std::unique_ptr<Link> >        LinkRefs;

class Solver{
public:
	struct Method{
		enum{
			GaussNewton,   ///< Gauss Newton method
			DDP,           ///< Differential Dynamic Programming
            Num
		};
	};
	
	struct Param{
		bool              verbose;
		int               method;
		std::vector<int>  numIter;          ///< number of minor iterations
		real_t            minStepSize;
		real_t            maxStepSize;
		real_t            cutoffStepSize;
		bool              hastyStepSize;
		real_t            regularization;
		real_t            stateRegularization;
        real_t            complRelaxation;  ///< relaxation parameter of barrier inequality constraints
		bool              fixInitialState;  ///< fix initial state in ddp
		bool              fixInitialInput;  ///< fix initial input in ddp
		bool              parallelize;
		
		Param();
	};

	struct Status{
		real_t  obj;        ///< value of objective function
		real_t  objDiff;    ///< change of value of objective function
		real_t  stepSize;   ///< step size
		int     iterCount;  ///< cumulative iteration count
		int     timePre;
		int     timeDir;
		int     timeStep;
		int     timeMod;
		int     timeTrans, timeTransRev;
		int     timeCost, timeCostGrad;
		int     timeBack;

		Status();
	};

	struct VariableInfo{
		int  num;
		int  numLocked;

		VariableInfo();
	};

	struct ConstraintInfo{
		int     num;
		int     numEnabled;
		int     numActive;
		real_t  error;

		ConstraintInfo();
	};

	struct Request{
		struct Type{
			enum{
				Enable,
				Lock,
				SetPriority,
				SetCorrection,
				SetConstraintWeight,
				SetVariableWeight,
			};
		};

		int    type;
		ID     mask;
		bool   enable;
		bool   lock;
		real_t rate;
		real_t lim;
		vec3_t weight;
	};

	struct SubState{
		int                 index;
		Variable*           var;
	};
	struct SubInput{
		int                 index;
		Variable*           var;
	};
	struct SubStateLink{
		SubState* x;
		Link*     link;
	};
	struct SubInputLink{
		SubInput* u;
		Link*     link;
	};
		
	struct SubTransition{
		Constraint*                 con;
		SubState*                   x1;
		std::vector<SubStateLink>   x0;
		std::vector<SubInputLink>   u;
	};
	struct SubCost{
		int                         index;
        Constraint*                 con;
		std::vector<SubStateLink>   x;
		std::vector<SubInputLink>   u;
		int  xbegin, xend;
		int  ubegin, uend;
	};

	struct State{
		int  dim;
		std::vector< std::unique_ptr<SubState> > substate;

		SubState*  Find(Variable* var);
	};
	struct Input{
		int  dim;
		std::vector< std::unique_ptr<SubInput> > subinput;

		SubInput*  Find(Variable* var);
	};
	struct Transition{
		std::vector< std::unique_ptr<SubTransition> > subtran;
	};
	struct Cost{
		int     dim;
		std::vector< std::unique_ptr<SubCost> >  subcost;
		Vector  y, b;
		Matrix  Ax;
		Matrix  Au;
	};
	
	Param            param;
	Status           status;
	std::vector<Request>  requests;
	
	VariableRefs         vars;			///< array of all variables
	Variables            vars_unlocked;
	ConstraintRefs       cons;			///< array of all constraints
	Constraints          cons_active;	///< array of active constraints
	
	LinkRefs        links;			///< array of links

	int         dimvar         ;
	int         dimvar_weighted;
	int         dimcon         ;
	Matrix      A, AtrA;
	Vector      b, b2;
	Vector      Atrb;
	Vector      dxvec;
	std::vector<int> pivot;

	std::vector<VariableInfo>    varInfoType;
	std::vector<ConstraintInfo>  conInfoType;		///< sum for each constraint category
	
	bool  ready;

	std::vector< std::unique_ptr<State> >       state;
	std::vector< std::unique_ptr<Input> >       input;
	std::vector< std::unique_ptr<Transition> >  transition;
	std::vector< std::unique_ptr<Cost> >        cost;

	int                    N;
	std::vector<real_t>    dt;
	std::vector<Vector>    dx;
	std::vector<Vector>    du;
	std::vector<Matrix>    fx;
	std::vector<Matrix>    fu;
	std::vector<Vector>    fcor;
    std::vector<real_t>    L;
	std::vector<Vector>    Lx;
	std::vector<Matrix>    Lxx;
	std::vector<Vector>    Lu;
	std::vector<Matrix>    Luu;
	std::vector<Matrix>    Lux;
	std::vector<real_t>    Q;
	std::vector<Vector>    Qx;
	std::vector<Vector>    Qu;
	std::vector<Matrix>    Qxx;
	std::vector<Matrix>    Quu;
	std::vector<Matrix>    Qux;
	std::vector<Matrix>    Quuinv;
	std::vector<Vector>    Quuinv_Qu;
	std::vector<real_t>    V, dV;
	std::vector<Vector>    Vx, dVx;
	std::vector<Matrix>    Vxx, dVxx;
	std::vector<Vector>    Vxx_fcor;
	std::vector<Matrix>    Vxx_fx;
	std::vector<Matrix>    Vxx_fu;
    std::vector<Vector>    Vx_plus_Vxx_fcor;
	std::vector<Matrix>    Quuinv_Qux;
	std::vector<Vector>    Qu_plus_Qux_dx;
	Matrix                 Vxxinv;

	Timer timer, timer2, timer3, timer4;
	
public:
	/// internal functions
	void    Prepare             ();
	real_t  CalcUpdatedObjective(real_t alpha);
	void    CalcEquation        ();
	void    InitDDP             ();
	void    ClearDDP            ();
	void    CalcTransitionDDP   ();
	void    CalcCostDDP         ();
	void    CalcCostGradientDDP ();
	void    BackwardDDP         ();
	void    ForwardDDP          (real_t alpha);
	void    CalcDirectionDDP    ();
	real_t  CalcObjectiveDDP    ();
	void    ModifyVariablesDDP  (real_t alpha);
	
	void AddVar   (Variable* var);      ///< add variable
	void DeleteVar(Variable* var);	    ///< delete variable
	void AddCon   (Constraint* con);    ///< add constraint
	void DeleteCon(Constraint* con);    ///< delte variable

	/// methods for DDP
	SubState*             AddStateVar            (Variable*   var, int k);  ///< register var as a sub-state at k
	SubInput*             AddInputVar            (Variable*   var, int k);  ///< register var as a sub-input at k
	SubTransition*        AddTransitionCon       (Constraint* con, int k);  ///< register con as a forward transition at k
	SubCost*              AddCostCon             (Constraint* con, int k);  ///< register con as a cost at k
	
public:
	/// virtual function that are to be overridden by derived classes

	/// 
	virtual void    ModifyVariables     (real_t alpha);

	/// evaluate objective function
	virtual real_t  CalcObjective();

	/// calculate update direction
	virtual void    CalcDirection();

	/// calculate step size
	virtual real_t  CalcStepSize();
	
public:
	/** @brief	enable or disable constraints
		@param	mask	constraint id mask
		@param	enable	enable or disable
	 */
	void Enable(ID mask, bool enable = true);

	/** @brief  lock or unlock variables
		@param	mask	variable id mask
		@param	lock	lock or unlock
	 */
	void Lock(ID mask, bool lock = true);

	/** @brief	set correction rate
	 */
	void SetCorrection(ID mask, real_t rate, real_t lim = FLT_MAX);

	/** @brief	set constraint weight
	    @param  mask    constraint id mask
		@param  weight  constraint weight

		Solver minimizes the weighted sum of constraint errors and variable changes.
		This function sets the weight of constraints specified by mask.
	 */
	void SetConstraintWeight(ID mask, real_t weight);
	void SetConstraintWeight(ID mask, vec3_t weight);
	
	/** @brief	set variable weight
	    @param  mask    variable id mask
		@param  weight  variable weight

		Solver minimizes the weighted sum of constraint errors and variable changes.
		This function sets the weight of variables specified by mask.
     */
	void SetVariableWeight  (ID mask, real_t weight);
	void SetVariableWeight  (ID mask, vec3_t weight);
	
	/** @brief calculate constraint error
		@param mask			constraint id mask
		@param sum_or_max	if true, sum of constraint errors is returned. otherwise the maximum is returned.
		
		for constraints with IDs that matches mask, sum of max of constraint error
		is calculated depending on sum_or_max.
	 */
	real_t CalcError(ID mask, bool sum_or_max);

	/// do initialization
	virtual void Init();

	/// one step update
	void Step();

	/// deletes all variables and constraints
	void Clear();

	void Reset();

	Solver();

};

}
