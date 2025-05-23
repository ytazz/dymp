﻿#pragma once

#include <dymp/model.h>
#include <dymp/variable.h>
#include <dymp/constraint.h>
#include <dymp/blas.h>
#include <dymp/timer.h>
#include <dymp/util.h>

namespace dymp{;

class  Wholebody;
struct WholebodyJointPosCon;
struct WholebodyJointVelCon;
struct WholebodyJointAccCon;
struct WholebodyCentroidPosConT;
struct WholebodyCentroidVelConT;
struct WholebodyCentroidPosConR;
struct WholebodyCentroidLCon;
struct WholebodyDesPosConT;
struct WholebodyDesPosConR;
struct WholebodyDesVelConT;
struct WholebodyDesVelConR;
struct WholebodyContactPosConT;
struct WholebodyContactPosConR;
struct WholebodyContactVelConT;
struct WholebodyContactVelConR;
struct WholebodyNormalForceCon;
struct WholebodyFrictionForceCon;
struct WholebodyMomentCon;

class WholebodyCallback;

struct WholebodyData{
	struct Link{
		vec3_t  pos_t;  // position (base link local)
		quat_t  pos_r;  // orientation
		vec3_t  vel_t;  // velocity
		vec3_t  vel_r;  // angular velocity
		vec3_t  acc_t;  // acceleration
		vec3_t  acc_r;  // angular acceleration
		vec3_t  force_t, force_t_par, force_t_child;
		vec3_t  force_r, force_r_par, force_r_child;
		mat3_t  I;

		Link();
	};

	struct End{
		vec3_t  pos_t, pos_t_abs;
		quat_t  pos_r, pos_r_abs;
		vec3_t  vel_t, vel_t_abs;
		vec3_t  vel_r, vel_r_abs;
		vec3_t  acc_t;
		vec3_t  acc_r;
		vec3_t  force_t;
		vec3_t  force_r;

		vec3_t  pos_t_weight, pos_r_weight;
		vec3_t  vel_t_weight, vel_r_weight;
		vec3_t  force_t_weight, force_r_weight;
		
		int     state;       ///< contact state
		real_t  mu;          ///< friction
		vec3_t  pos_tc;      ///< contact point (end local)
		vec3_t  normal;      ///< contact normal (global)
		vec3_t  cop_min;
		vec3_t  cop_max;

		End();
	};

    struct Joint{
        real_t q, qd, qdd, qddd, tau;
	    real_t q_weight, qd_weight, qdd_weight, qddd_weight;
        real_t q_range_weight, qd_range_weight, qdd_range_weight;
    };

	struct Centroid{
		vec3_t pos_t;   ///< com position
		quat_t pos_r;
		vec3_t vel_t;   ///< com velocity
		vec3_t vel_r;
		vec3_t acc_t;
		vec3_t acc_r;
		vec3_t pos_t_weight;
		vec3_t pos_r_weight;
		vec3_t vel_t_weight;
        vec3_t L_weight;
		vec3_t acc_t_weight;
        vec3_t Ld_weight;
		vec3_t L_local, Ld_local, L_abs;                     ///< momentum (local) and its derivative
		mat3_t I_local, Id_local, I_abs, Id_abs, I_abs_inv;  ///< inertia matrix around com and its inverse

		Centroid();
	};

	Centroid           centroid;
	std::vector<End>   ends;
	std::vector<Link>  links;
    std::vector<Joint> joints;
	
	Matrix          Jcom;
	std::vector<Matrix>  Jfk;
	
	void Init        (Wholebody* wb);
	void InitJacobian(Wholebody* wb);
	void CopyVars    (WholebodyData& d);
};

/*
 *  whold-body kinematics and (semi)dynamics
 */
class WholebodyKey : public Keypoint {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Wholebody*  wb;
	
	struct End{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		V3Var*  var_force_t;  ///< force (contact frame)
		V3Var*  var_force_r;  ///< moment

		WholebodyDesPosConT*  con_des_pos_t  ;   ///< desired position (global)
		WholebodyDesPosConR*  con_des_pos_r  ;   ///< desired orientation
		WholebodyDesVelConT*  con_des_vel_t  ;   ///< desired velocity
		WholebodyDesVelConR*  con_des_vel_r  ;   ///< desired angular velocity
		FixConV3*             con_des_force_t;   ///< desired force (contact frame)
		FixConV3*             con_des_force_r;   ///< desired moment
		
		WholebodyContactPosConT*    con_contact_pos_t;
		WholebodyContactPosConR*    con_contact_pos_r;
		WholebodyContactVelConT*    con_contact_vel_t;
		WholebodyContactVelConR*    con_contact_vel_r;
		WholebodyNormalForceCon*    con_force_normal;
		WholebodyFrictionForceCon*  con_force_friction[2][2];
		WholebodyMomentCon*         con_moment[3][2];
	};

	struct Centroid{		
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		V3Var*  var_pos_t;
		QVar*   var_pos_r;
		V3Var*  var_vel_t;
        V3Var*  var_L;
		V3Var*  var_acc_t;
        V3Var*  var_Ld;

		WholebodyCentroidPosConT*  con_pos_t;
		WholebodyCentroidPosConR*  con_pos_r;
		WholebodyCentroidVelConT*  con_vel_t;
        WholebodyCentroidLCon*     con_L;

		FixConV3*  con_des_pos_t;
		FixConQ*   con_des_pos_r;
		FixConV3*  con_des_vel_t;
        FixConV3*  con_des_L;
		FixConV3*  con_des_acc_t;
        FixConV3*  con_des_Ld;

	};

	struct Joint{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		SVar*  var_q;
		SVar*  var_qd;
		SVar*  var_qdd;
		SVar*  var_qddd;
		
		WholebodyJointPosCon*  con_q;
		WholebodyJointVelCon*  con_qd;
		WholebodyJointAccCon*  con_qdd;

		FixConS*  con_des_q;
		FixConS*  con_des_qd;
		FixConS*  con_des_qdd;
		FixConS*  con_des_qddd;

		RangeConS* con_range_q;
		RangeConS* con_range_qd;
		RangeConS* con_range_qdd;
	};

	WholebodyData          data;
	WholebodyData          data_des;

	Centroid            centroid;
	std::vector<End>    ends;
	std::vector<Joint>  joints;

	// working variables
	//quat_t          q0;
	//mat3_t          R0;
	std::vector<vec3_t>  re, fe, me;
	std::vector<mat3_t>  rec;
	vec3_t          fsum, msum;
	Matrix          J_L_q, J_L_qd;
	Matrix          J_Ld_q, J_Ld_qdd;
	Matrix          mj_pjc, mj_vjc, mj_ajc, Ij;
	
public:	
    virtual void AddVar(Solver* solver);
	virtual void AddCon(Solver* solver);
	virtual void Prepare     ();
	virtual void PrepareStep ();
	virtual void Finish      ();
	
	WholebodyKey();
};

class Wholebody : public Model{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	struct ContactState{
		enum{
			Free,
			Surface,
			Line,
			Point,
		};
	};
	struct InputMode{
		enum{
			Velocity,
			Acceleration,
			Jerk,
		};
	};

	struct Param {
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		real_t  totalMass;  ///< total mass of wholebody
		vec3_t  nominalInertia;
		real_t  gravity;
		real_t  dt;         ///< time resolution. used for scaling only
		bool    useLd;
		//bool    useJerk;
		int     inputMode;
		
		Param();
	};

	struct Scale{
		real_t  l;
		real_t  t, tinv;   //< time scaling
		real_t  pt;  //< position scaling
		real_t  pr;
		real_t  vt;  //< velocity scaling
		real_t  vr;  //< angular velocity scaling
		real_t  at;
		real_t  ar;
		real_t  jt;
		real_t  jr;
		real_t  ft;  //< force scaling
		real_t  fr;  //< moment scaling
		real_t  L;   //< momentum scaling
	};

	struct Joint{
		real_t		rotor_inertia;
        vec2_t      pos_range;
        vec2_t      vel_range;
        vec2_t      acc_range;

		Joint(real_t Ir = 0.0);
	};

	struct Link{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		real_t       mass;         ///< mass of link
		real_t       mass_ratio;
		mat3_t       inertia, inertia_with_rotor;
		vec3_t       center;
		int          iparent;      ///< parent link index
		std::vector<int>  ichildren;    ///< child link indices
		int          ijoint;       ///< joint index
		int          iend;
		vec3_t       trn;          ///< translation from parent
		vec3_t       axis;         ///< joint axis
	
		Link(real_t _mass = 0.0, const vec3_t& _inertia = zero3, const vec3_t& _center = zero3, int _iend = -1, int _iparent = -1, int _ijoint = -1, const vec3_t& _trn = zero3, const vec3_t& _axis = zero3);
	};

	struct End{
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		int    ilink;    ///< link index
		vec3_t offset;
		bool   enableTranslation;
		bool   enableRotation;
		bool   enableForce;
		bool   enableMoment;
        bool   enableTerminalCost;

		End(int _ilink = 0.0, const vec3_t& _offset = zero3, bool _enable_trn = true, bool _enable_rot = true, bool _enable_force = true, bool _enable_moment = true, bool _enable_terminal = true);
	};
	
   	struct Snapshot{
		struct Link{
			vec3_t  pos_t;
			quat_t  pos_r;
			vec3_t  vel_t;
			vec3_t  vel_r;
			vec3_t  force_t;
			vec3_t  force_r;
		};
	
		real_t  t;
		vec3_t  pos_t;
		vec3_t  vel_t;
		quat_t  pos_r;
		vec3_t  vel_r;
		std::vector<Link>  links;
		
		Snapshot();
	};

    Param	            param;
	Scale               scale;
    std::vector<Link>   links;
	std::vector<Joint>  joints;
	std::vector<End>    ends;
	WholebodyCallback*  callback;
	
	Snapshot                 snapshot;
	std::vector<Snapshot>    trajectory;
	bool                     trajReady;

	Timer timer;
	Timer timer2;

	virtual Keypoint*	CreateKeypoint() { return new WholebodyKey(); }
	virtual void		Init   ();
	virtual void		Prepare();
	virtual void		PrepareStep();
	virtual void        Finish ();
	virtual void        CreateSnapshot(real_t time);
	
	void SetScaling();
	void Reset(bool reset_all);
	void Shift(real_t offset);
	void Setup();
	void CalcFK                (WholebodyData& d);
	void CalcPosition          (WholebodyData& d);
	void CalcJacobian          (WholebodyData& d);
	void CalcVelocity          (WholebodyData& d);
	void CalcAcceleration      (WholebodyData& d);
	void CalcComAcceleration   (WholebodyData& d);
	void CalcBaseAngularVelocity    (WholebodyData& d);
	void CalcBaseAngularAcceleration(WholebodyData& d);
	void CalcInertia                (WholebodyData& d);
	void CalcInertiaDerivative      (WholebodyData& d);
	void CalcLocalMomentum          (WholebodyData& d);
	void CalcLocalMomentumDerivative(WholebodyData& d);
	void CalcAbsoluteMomentum       (WholebodyData& d);
	void CalcForce                  (WholebodyData& d);
	
	void ComState    (real_t t, vec3_t& pos, vec3_t& vel   );
	void BaseState   (real_t t, quat_t& ori, vec3_t& angvel);
	void LinkPose    (real_t t, int i, vec3_t& pos  , quat_t& ori   );
    void LinkVelocity(real_t t, int i, vec3_t& vel  , vec3_t& angvel);
    void LinkForce   (real_t t, int i, vec3_t& force, vec3_t& moment);
    
	void CreateSnapshot(real_t t, Snapshot& s);
	void CalcTrajectory();
	
public:
	         Wholebody(World* g, string n);
	virtual ~Wholebody();
};

class WholebodyCallback{
public:
	virtual void   GetInitialState(WholebodyData& d) = 0;
	virtual void   GetDesiredState(int k, real_t t, WholebodyData& d) = 0;
};

struct WholebodyCon : Constraint {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	WholebodyKey*  obj[2];

	WholebodyCon(Solver* solver, int _dim, int _tag, string _name, WholebodyKey* _obj, real_t _scale);
};

struct WholebodyJointPosCon : WholebodyCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int    ijoint;
	real_t q0, qd0, qdd0, qddd0, q1, q_rhs;
	real_t h, h2, h3;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyJointPosCon(Solver* solver, string _name, WholebodyKey* _obj, int _ijoint, real_t _scale);
};

struct WholebodyJointVelCon : WholebodyCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int    ijoint;
	real_t qd0, qdd0, qddd0, qd1, qd_rhs;
	real_t h, h2;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyJointVelCon(Solver* solver, string _name, WholebodyKey* _obj, int _ijoint, real_t _scale);
};

struct WholebodyJointAccCon : WholebodyCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int    ijoint;
	real_t qdd0, qddd0, qdd1, qdd_rhs;
	real_t h;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyJointAccCon(Solver* solver, string _name, WholebodyKey* _obj, int _ijoint, real_t _scale);
};

struct WholebodyCentroidPosConT : WholebodyCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	real_t h, h2, m;
	vec3_t g;
	vec3_t pc0, vc0, ac0, pc1, pc_rhs;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyCentroidPosConT(Solver* solver, string _name, WholebodyKey* _obj, real_t _scale);
};

struct WholebodyCentroidVelConT : WholebodyCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	real_t h, m;
	vec3_t g;
	vec3_t vc0, ac0, vc1, vc_rhs;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyCentroidVelConT(Solver* solver, string _name, WholebodyKey* _obj, real_t _scale);
};

struct WholebodyCentroidPosConR : WholebodyCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	quat_t q0, q1, q_rhs, q_omega;
	vec3_t w0, u0, L, Ld, omega;
	mat3_t Id, Iinv;
	real_t h, h2;
	mat3_t R0, R_omega, A_omega;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyCentroidPosConR(Solver* solver, string _name, WholebodyKey* _obj, real_t _scale);
};
struct WholebodyCentroidLCon : WholebodyCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	vec3_t L0, L1, Ld0, L_rhs;
	real_t h;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyCentroidLCon(Solver* solver, string _name, WholebodyKey* _obj, real_t _scale);
};

struct WholebodyDesPosConT : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	WholebodyKey*  obj;
	int    iend;
	vec3_t desired;
	vec3_t pc, pe, pe_abs, oe, pi, ci, r;
	quat_t q0, qi;
	mat3_t R0;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyDesPosConT(Solver* solver, string _name, WholebodyKey* _obj, int _iend, real_t _scale);
};

struct WholebodyDesPosConR : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	WholebodyKey*  obj;
	int    iend;
	quat_t desired;
	quat_t q0, qe, qe_abs;
	mat3_t R0;

	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyDesPosConR(Solver* solver, string _name, WholebodyKey* _obj, int _iend, real_t _scale);
};

struct WholebodyDesVelConT : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	WholebodyKey*  obj;
	int    iend;
	vec3_t desired;
	vec3_t vc, w0, ve, ve_abs, oe, pi, ci, r, pe;
	quat_t q0, qi;
	mat3_t R0, Iinv;

	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyDesVelConT(Solver* solver, string _name, WholebodyKey* _obj, int _iend, real_t _scale);
};

struct WholebodyDesVelConR : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	WholebodyKey*  obj;
	int    iend;
	vec3_t desired;
	vec3_t w0, we, we_abs;
	quat_t q0;
	mat3_t R0, Iinv;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
		
	WholebodyDesVelConR(Solver* solver, string _name, WholebodyKey* _obj, int _iend, real_t _scale);
};

struct WholebodyContactCon : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	WholebodyKey*  obj;
	int    iend;
	vec3_t po, r;
	quat_t qo;
	mat3_t Ro;

	void Prepare();

	WholebodyContactCon(Solver* solver, string _name, WholebodyKey* _obj, int _tag, int _iend, real_t _scale);
};

struct WholebodyContactPosConT : WholebodyContactCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	vec3_t pc, pe;
	quat_t q0, qe;
	mat3_t R0;

	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
	
	WholebodyContactPosConT(Solver* solver, string _name, WholebodyKey* _obj, int _iend, real_t _scale);
};

struct WholebodyContactPosConR : WholebodyContactCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	quat_t q0, qe;
	mat3_t R0;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
	
	WholebodyContactPosConR(Solver* solver, string _name, WholebodyKey* _obj, int _iend, real_t _scale);
};

struct WholebodyContactVelConT : WholebodyContactCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	vec3_t vc, w0, pe, ve, we;
	quat_t q0, qe;
	mat3_t R0, Iinv;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
	
	WholebodyContactVelConT(Solver* solver, string _name, WholebodyKey* _obj, int _iend, real_t _scale);
};

struct WholebodyContactVelConR : WholebodyContactCon{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	quat_t q0, qe;
	mat3_t R0, Iinv;
	vec3_t w0, we;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
	
	WholebodyContactVelConR(Solver* solver, string _name, WholebodyKey* _obj, int _iend, real_t _scale);
};

struct WholebodyNormalForceCon : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	WholebodyKey*  obj;
	int    iend;
	quat_t qi;
	vec3_t nz, f;
	real_t fz;

	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
	
	WholebodyNormalForceCon(Solver* solver, string _name, WholebodyKey* _obj, int _iend, real_t _scale);
};

struct WholebodyFrictionForceCon : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	WholebodyKey*  obj;
	int    iend;
	int    dir;   //< x or y
	int    side;  //< upper or lower bound
	quat_t qi;
	vec3_t nx, ny, nz, f, df;
	real_t fz, ft, mu;
	
	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
	
	WholebodyFrictionForceCon(Solver* solver, string _name, WholebodyKey* _obj, int _iend, int _dir, int _side, real_t _scale);
};

struct WholebodyMomentCon : Constraint{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	WholebodyKey*  obj;
	int    iend;
	int    dir;   //< x or y
	int    side;  //< upper or lower bound
	quat_t qi;
	vec3_t nx, ny, nz, f, eta, df, deta;
	real_t fz, etax, etay, etaz;
	vec3_t cmin, cmax;

	void Prepare();

	virtual void  CalcCoef();
	virtual void  CalcDeviation();
	
	WholebodyMomentCon(Solver* solver, string _name, WholebodyKey* _obj, int _iend, int _dir, int _side, real_t _scale);
};

}
