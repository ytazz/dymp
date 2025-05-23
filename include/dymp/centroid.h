﻿#pragma once

#include <dymp/model.h>
#include <dymp/variable.h>
#include <dymp/constraint.h>

namespace dymp{;

class  Centroid;
struct CentroidPosConT;
struct CentroidPosConR;
struct CentroidVelConT;
struct CentroidLCon;
struct CentroidTimeCon;
struct CentroidEndPosConT;
struct CentroidEndPosConR;
struct CentroidDurationRangeCon;
struct CentroidEndPosRangeCon;
//struct CentroidEndContactPosConT;
//struct CentroidEndContactPosConR;
//struct CentroidEndContactVelConT;
//struct CentroidEndContactVelConR;
struct CentroidEndFrictionCon;
struct CentroidEndMomentCon;

class CentroidCallback;

struct CentroidData{
	struct End{
		/// end pose in contact surface coordinate
		vec2_t  pos_t;         //< end position
		real_t  pos_r;         //< end orientation
		vec2_t  vel_t;         //< end velocity
		real_t  vel_r;         //< end angular velocity

		/// contact wrench (stiffness-based parametrization)
		real_t  stiff;         //< stiffness
		vec2_t  cop;           //< cop of stiffness-based parametrization
		vec2_t  cmp;           //< cmp of stiffness-based parametrization
		real_t  torsion;       //< torsion (z component of moment) of stiffness-based parametrization
		
		/// contact wrench (direct parametrization)
		vec3_t  force;         //< force value of direct parametrization
		vec3_t  moment;        //< moment value of direct parametrization
		
		/// auxiliary contact info
		vec3_t  pos_tc;        //< contact point (end-local)
		int     iface;         //< contact face index
		bool    contact;       //< contact state
		real_t  mu;
		vec3_t  cop_min;       //< cop range
		vec3_t  cop_max;
		real_t  roll;
		real_t  rolld;
		real_t  tilt;          //< tilt angle
		real_t  tiltd;         //< tilt angular velocity
		real_t  alt;
		real_t  altd;

		/// decoded 3D info
		vec3_t  pos_t_abs;
		quat_t  pos_r_abs;
		vec3_t  vel_t_abs;
		vec3_t  vel_r_abs;
		vec3_t  force_t;
		vec3_t  force_r;
		
		/// auxiliary info for interpolation
		real_t  t_lift, t_land;
		int     iface_lift, iface_land;
		vec3_t  pos_t_lift, pos_t_land;
		quat_t  pos_r_lift, pos_r_land;
		vec3_t  pos_t_lift_local, pos_t_land_local;
		quat_t  pos_r_lift_local, pos_r_land_local;
		vec3_t  vel_t_ave, vel_r_ave;
		//vec2_t  vel_t_ave;
		//real_t  vel_r_ave;
		
		/// weights
		vec2_t  pos_t_weight;
		real_t  pos_r_weight;
		vec2_t  vel_t_weight;
		real_t  vel_r_weight;
		real_t  stiff_weight;
		vec2_t  cop_weight;
		vec2_t  cmp_weight;
		real_t  torsion_weight;
		vec3_t  force_weight;
		vec3_t  moment_weight;
		
		End();
	};

	vec3_t pos_t;
	quat_t pos_r;
	vec3_t vel_t;
	vec3_t vel_r;
	vec3_t L;
	vec3_t acc_t;
	vec3_t acc_r;
	real_t time;
	real_t duration;
	
	vec3_t pos_t_weight;
	vec3_t pos_r_weight;
	vec3_t vel_t_weight;
	vec3_t L_weight;
	real_t time_weight;
	real_t duration_weight;

	// inertia
	std::vector<mat3_t> I, Iinv;
	std::vector<vec3_t> Llocal;

	// auxiliary data
	real_t lbar;
	vec3_t pbar, rbar, etabar;
	vec3_t fsum, etasum;

	vec3_t              p_rhs;
	std::vector<vec3_t> v_rhs;	
	std::vector<vec3_t> w_rhs;
	std::vector<vec3_t> L_rhs;
	std::vector<quat_t> q_rhs;
	
	std::vector<End>   ends;
		
	void Init(Centroid* cen);
	void CopyVars(CentroidData& d);

	CentroidData();
};

/**
	centroidal dynamics model
*/
class CentroidKey : public Keypoint {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Centroid*     cen;
	CentroidKey*  endNext;

	CentroidData  data, data_des;

	int   iphase;   //< phase index
	int   idiv;     //< subdivision index

	V3Var*  var_pos_t;    ///< position        
	QVar*   var_pos_r;    ///< orientation
	V3Var*  var_vel_t;    ///< velocity        
	V3Var*  var_L;
	SVar*   var_time;
	SVar*   var_duration;  //< duration of this contact phase
	
	CentroidPosConT*           con_pos_t;
	CentroidPosConR*           con_pos_r;
	CentroidVelConT*           con_vel_t;
	CentroidLCon*              con_L;
	CentroidTimeCon*           con_time ;
	CentroidDurationRangeCon*  con_duration_range[2];
    FixConV3*                  con_des_pos_t;
    FixConQ*                   con_des_pos_r;
    FixConV3*                  con_des_vel_t;
	FixConV3*                  con_des_L;
	FixConS*                   con_des_time;
	FixConS*                   con_des_duration;

	struct End{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		CentroidKey*  key;
		
		SVar*  var_pos_t[2];  //< end effector position (in global coordinate)
		SVar*  var_pos_r;
		SVar*  var_vel_t[2];  //< end effector velocity (in global coordinate)
		SVar*  var_vel_r;
		SVar*  var_stiff;     //< contact stiffness
		SVar*  var_cop[2];    //< contact cop
		SVar*  var_cmp[2];    //< contact cmp
		SVar*  var_torsion;   //< contact torsion
		V3Var* var_force;     //< contact force
		V3Var* var_moment;    //< contact moment

		CentroidEndPosConT*             con_pos_t[2];
		CentroidEndPosConR*             con_pos_r;
		CentroidEndPosRangeCon*         con_pos_range[3][2];
		RangeConS*                      con_stiff_range;
		FixConS*                        con_des_pos_t[2];
		FixConS*                        con_des_pos_r;
		FixConS*                        con_des_vel_t[2];
		FixConS*                        con_des_vel_r;
		FixConS*                        con_des_stiff;
		FixConS*                        con_des_cop[2];
		FixConS*                        con_des_cmp[2];
		FixConS*                        con_des_torsion;
		FixConV3*                       con_des_force;
		FixConV3*                       con_des_moment;

		CentroidEndFrictionCon*         con_friction[2][2];
		CentroidEndMomentCon*           con_moment[3][2];		
	};

	vector<End>    ends;

	vec3_t  g;
	real_t  m;
	real_t  dtau;
	vector<real_t> t, t2, t3, C, S;
	vector<real_t> C_tau, S_tau, C_lbar, S_lbar;
	vector<real_t> li, li2;
	vector<vec3_t> pi, ci, ri, ti, fi, etai;
	
	vec3_t psum, rsum;
	real_t l2sum;
	mat3_t pbar_cross, rbar_cross;
	
	vector<real_t> lbar_li;
	vec3_t         pbar_lbar;
	vector<vec3_t> pbar_li;
	vector<real_t> pbar_pi, pbar_ci;
	vec3_t         rbar_lbar;
	vector<vec3_t> rbar_li;
	vector<real_t> rbar_ri;
	vec3_t         etabar_lbar;
	mat3_t         etabar_pbar, etabar_rbar;
	vector<real_t> etabar_ti;
	vector<vec3_t> etabar_li;
	vector<mat3_t> etabar_pi, etabar_ci, etabar_ri;
	vector<mat3_t> R_omega;
	
	vec3_t p_C, p_S;
	real_t p_p, p_v, p_pbar, p_rbar;
	vec3_t p_lbar, p_tau;
	vector<vec3_t> p_li;              //< nend
	vector<real_t> p_pi, p_ci, p_ri, p_fi;  //< nend
	
	vec3_t v_C, v_S;
	vector<real_t> v_p, v_v, v_pbar, v_rbar;     //< ndiv array
	vector<vec3_t> v_lbar, v_tau;                //< ndiv array
	vector< vector<vec3_t> >  v_li;              //< nend x ndiv array
	vector< vector<real_t> >  v_pi, v_ci, v_ri, v_fi;  //< nend x ndiv array
	
	real_t L_L;
	vector<mat3_t> L_v1, L_p, L_v;                             //< ndiv
	vector<vec3_t> L_tau;                                //< ndiv
	vector<mat3_t> L_rbar;                               //< ndiv
	vector<real_t> L_etabar;                             //< ndiv
	vector< vector<mat3_t> >  L_pi, L_ci, L_ri, L_fi;    //< nend x ndiv
	vector< vector<real_t> >  L_etai, L_ti;
	vector< vector<vec3_t> >  L_li;                      //< nend x ndiv

	mat3_t q_q, q_L;
	mat3_t q_p, q_v;
	vec3_t q_tau;
	vector<mat3_t> q_L1;                             //< ndiv
	vector<mat3_t> q_pi, q_ci, q_ri, q_ti, q_fi, q_etai;   //< nend
	vector<vec3_t> q_li;                             //< nend
	
public:
	void Prepare2();

    virtual void AddVar(Solver* solver);
	virtual void AddCon(Solver* solver);
	virtual void Prepare();
	virtual void PrepareStep();
	virtual void Finish ();
	virtual void Draw(render::Canvas* canvas, render::Config* conf);

	CentroidKey();
};

class Centroid : public Model{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	/*
	struct ContactState{
		enum{
			Free,
			Surface,
			Line,
			Point,
		};
	};
	*/
	struct EndInterpolation{
		enum{
			Polynomial,
			CycloidLocal,
			CycloidGlobal,
		};
	};
	struct EndWrenchParametrization{
		enum{
			Direct,
			Stiffness,
		};
	};
	struct Param {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		real_t	m;  //< mass
		mat3_t  I;  //< inertia
		mat3_t  Iinv;
		real_t	g;  //< gravity
		//real_t  mu;

        real_t  durationMin;
		real_t  durationMax;
        
		real_t  swingSlope;
        real_t  swingHeight;
        real_t  swingLiftMargin;
		real_t  swingLandMargin;
		real_t  swingTurnMargin;

        bool    enableRotation;    ///< enable rotational dynamics
		int     rotationResolution;
		int     endInterpolation;
		int     endWrenchParametrization;

		real_t  complWeight;

		//real_t  contactMargin;
		
        Param();
	};
	
	struct Scale{
		real_t  l;
		real_t  t, tinv;   //< time scaling
		real_t  pt, pt2;  //< position scaling
		real_t  pr;
		real_t  vt;  //< velocity scaling
		real_t  vr;  //< angular velocity scaling
		real_t  at;
		real_t  ar;
		real_t  ft;  //< force scaling
		real_t  fr;  //< moment scaling
		real_t  L;   //< momentum scaling
	};

   	struct End{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		vec3_t  basePos;
		vec3_t  posMin;
		vec3_t  posMax;
		vec3_t  copMin;
        vec3_t  copMax;
        real_t  stiffnessMax;
        bool    lockOri;
		
        End();
	};

	struct Phase{
		vector<int>  iface;    ///< contact face index of each end. -1 means no contact
		int  ndiv;             ///< number of subdivisions
	};

	struct Waypoint {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		struct End{
			EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			struct Value{
				EIGEN_MAKE_ALIGNED_OPERATOR_NEW
				vec2_t  pos_t;  //< specified in local coordinate
				real_t  pos_r;
				vec2_t  vel_t;
				real_t  vel_r;
				int     iface;

				Value();
				Value(const vec2_t& _p, real_t _q, const vec2_t& _v, real_t _w, int _iface);
			};
			struct Weight{
				EIGEN_MAKE_ALIGNED_OPERATOR_NEW
				vec2_t  pos_t ;
				real_t  pos_r ;
				vec2_t  vel_t ;
				real_t  vel_r ;
				real_t  stiff ;
				vec2_t  cop;
				vec2_t  cmp;
				real_t  torsion;
				vec3_t  force;
				vec3_t  moment;

				Weight();
				Weight(const vec2_t& _p, real_t _q, const vec2_t& _v, real_t _w, real_t _l, const vec2_t& _c, const vec2_t& _r, real_t _t, const vec3_t& _f, const vec3_t& _m);
			};

			Value   value;
			Weight  weight;
			
            End();
		};

		struct Value{
			EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			real_t  time;
			real_t  duration;
			vec3_t  pos_t;
			vec3_t  pos_r;
			vec3_t  vel_t;
			vec3_t  vel_r;

			Value();
			Value(real_t _t, real_t _tau, const vec3_t& _p, const vec3_t& _q, const vec3_t& _v, const vec3_t& _w);
		};
		struct Weight{
			EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			real_t  time;
			real_t  duration;
			vec3_t  pos_t;
			vec3_t  pos_r;
			vec3_t  vel_t;
			//vec3_t  vel_r;
			vec3_t  L;

			Weight();
			Weight(real_t _t, real_t _tau, const vec3_t& _p, const vec3_t& _q, const vec3_t& _v, const vec3_t& _L);//_w);
		};

		Value   value;
		Weight  weight;
		
		vector<End>  ends;

		Waypoint();
	};
	
    struct Face{
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		vec3_t               pos_t;
		quat_t               pos_r;
		vec3_t               normal;
		vec2_t               slope;
		real_t               elevation;
		std::vector<vec3_t>  vertices;

		void Init();
        
        Face();
    };

    Param	            param;
	Scale               scale;
    vector<End>         ends;
    vector<Face>        faces;
	vector<Phase>       phases;
	vector<Waypoint>    waypoints;
    
	CentroidData         snapshot;
	vector<CentroidData> trajectory;
	bool                 trajReady;

	CentroidCallback*   callback;

	int  NumSteps();
	void Reset(bool reset_first, bool reset_middle, bool reset_last);
	void Shift();
	void Setup();
	void CalcComAcceleration (CentroidData& d);
	void CalcBaseAcceleration(CentroidData& d);
	void CalcEndState        (int iend, CentroidData& d, bool _2d_to_3d);
	void CalcEndWrench       (int iend, CentroidData& d);
	void ResetEndWrench      (CentroidData& d);
	int  FindFace            (const vec3_t& p, const vec3_t& nz);
	int  CreateFace          (const vec3_t& p, const vec3_t& nz);
	
	virtual Keypoint*	CreateKeypoint() { return new CentroidKey(); }
	virtual void		Init   ();
	virtual void		Prepare();
	virtual void		PrepareStep();
	virtual void        Finish ();
	virtual void        CreateSnapshot(real_t time);
    virtual void        DrawSnapshot(render::Canvas* canvas, render::Config* conf);
	virtual void        Draw        (render::Canvas* canvas, render::Config* conf);
	
    void SetScaling     ();
	void CalcState      (real_t t, CentroidData& d);
	void CalcState      (real_t t, const vector<CentroidData>& d_array, CentroidData& d);
	void CalcState      (real_t t, const CentroidData& d0, const CentroidData& d1, CentroidData& d);
	//void CalcEndState   (real_t t, int iend, bool& contact, int& iface, vec2_t& pos_lift, real_t& ori_lift, vec2_t& pos_land, real_t& ori_land);
	void CalcTrajectory ();
	
public:
	         Centroid(World* w, string n);
	virtual ~Centroid();
};

class CentroidCallback{
public:
	virtual void   GetInitialState(CentroidData& d) = 0;
	virtual void   GetDesiredState(int k, real_t t, CentroidData& d) = 0;
};

struct CentroidCon : Constraint {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	CentroidKey*  obj[2];
	
	void Prepare();

	CentroidCon(Solver* solver, int _dim, int _tag, string _name, CentroidKey* _obj, real_t _scale);
};

struct CentroidComCon : CentroidCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	void Prepare();

	CentroidComCon(Solver* solver, int _dim, int _tag, string _name, CentroidKey* _obj, real_t _scale);
};

struct CentroidPosConT : CentroidComCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	void Prepare();
	
	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	CentroidPosConT(Solver* solver, string _name, CentroidKey* _obj, real_t _scale);
};

struct CentroidVelConT : CentroidComCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	CentroidVelConT(Solver* solver, string _name, CentroidKey* _obj, real_t _scale);
};

struct CentroidPosConR : CentroidCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	CentroidPosConR(Solver* solver, string _name, CentroidKey* _obj, real_t _scale);
};
/*
struct CentroidVelConR : CentroidCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	CentroidVelConR(Solver* solver, string _name, CentroidKey* _obj, real_t _scale);
};
*/
struct CentroidLCon : CentroidCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	CentroidLCon(Solver* solver, string _name, CentroidKey* _obj, real_t _scale);
};

struct CentroidTimeCon : CentroidCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	real_t t, t_lhs, t_rhs;

	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	CentroidTimeCon(Solver* solver, string _name, CentroidKey* _obj, real_t _scale);
};

struct CentroidEndPosConT : CentroidCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int    iend;
	int    dir;
	real_t pe, ve, pe_lhs, pe_rhs;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	CentroidEndPosConT(Solver* solver, string _name, CentroidKey* _obj, int _iend, int _dir, real_t _scale);
};

struct CentroidEndPosConR : CentroidCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int    iend;
	real_t qe, we, qe_lhs, qe_rhs;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	CentroidEndPosConR(Solver* solver, string _name, CentroidKey* _obj, int _iend, real_t _scale);
};

struct CentroidDurationRangeCon : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	CentroidKey* obj;
	real_t       dir;
	real_t       bound;
	
	virtual void  CalcCoef();
	virtual void  CalcDeviation();

	CentroidDurationRangeCon(Solver* solver, string _name, CentroidKey* _obj, real_t _dir, real_t _scale);
};

struct CentroidEndPosRangeCon : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	CentroidKey* obj;
	int          iend;
	int          dir;
	int          side;
	vec3_t       p, pbase, pe, pf, eta, eta_abs;
	vec2_t       pe_local;
	quat_t       q, qf;
	mat3_t       Rf;
	real_t       bound;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();

	CentroidEndPosRangeCon(Solver* solver, string _name, CentroidKey* _obj, int _iend, int _dir, int _side, real_t _scale);
};

struct CentroidEndFrictionCon : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	CentroidKey*  obj;
	int           iend;
	int           dir;   //< x or y
	int           side;  //< upper or lower bound
    quat_t        qf;
	vec3_t        nx, ny, nz, f, df;
	real_t        fz, ft, mu;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();

	CentroidEndFrictionCon(Solver* solver, string _name, CentroidKey* _obj, int _iend, int _dir, int _side, real_t _scale);
};

struct CentroidEndMomentCon : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	CentroidKey*  obj;
	int           iend;
    int           dir;   //< x or y
	int           side;  //< upper or lower bound
	vec2_t        nx, ny, c, dc;
	vec3_t        nz, p, pf, f, eta, dp, df, deta;
	real_t        t, dt, cmin, cmax, bound;

	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();

	CentroidEndMomentCon(Solver* solver, string _name, CentroidKey* _obj, int _iend, int _dir, int _side, real_t _scale);
};

}
