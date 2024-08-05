#pragma once

#include <types.h>
#include <blas.h>

#include <limits>
using namespace std;

using namespace Eigen;

namespace dimp3{;

const real_t eps     = 1.0e-10;
const real_t _pi     = 3.1415926535;
const real_t inf     = numeric_limits<real_t>::max();
const int    infi    = numeric_limits<int>::max();
const vec2_t inf2    = vec2_t(inf, inf);
const vec3_t inf3    = vec3_t(inf, inf, inf);
const vec2_t zero2(0.0, 0.0);
const vec3_t zero3(0.0, 0.0, 0.0);
const vec2_t one2(1.0, 1.0);
const vec3_t one3(1.0, 1.0, 1.0);
const vec3_t ex  (1.0, 0.0, 0.0);
const vec3_t ey  (0.0, 1.0, 0.0);
const vec3_t ez  (0.0, 0.0, 1.0);

inline real_t deg_to_rad(real_t deg){
    return (_pi/180.0)*deg;
}

inline real_t rad_to_deg(real_t rad){
    return (180.0/_pi)*rad;
}

inline real_t square(real_t x){
	return x*x;
}

inline mat3_t vvtrmat(vec3_t c, vec3_t r){
	mat3_t m;
	for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++)
		m(i,j) = c[i]*r[j];
	return m;
}

inline mat3_t cross_mat(vec3_t r){
    mat3_t m;
    m <<  0   , -r(2),  r(1),
          r(2),  0   , -r(0),
         -r(1),  r(0),  0;
    return m;
}

inline void cross_mat(vec3_t r, real_t k, Matrix& m){
	m(0,0) =  0.0    ; m(0,1) = -k*r.z(); m(0,2) =  k*r.y();
	m(1,0) =  k*r.z(); m(1,1) =  0.0    ; m(1,2) = -k*r.x();
	m(2,0) = -k*r.y(); m(2,1) =  k*r.x(); m(2,2) =  0.0    ;
}

inline quat_t rot_quat(vec3_t w){
    real_t wnorm = w.norm();
    if(wnorm < eps)
        return quat_t(1.0, 0.0, 0.0, 0.0);

    return quat_t(AngleAxisd(wnorm, w/wnorm));
}

inline vec3_t quat_error(quat_t q0, quat_t q1){
    AngleAxisd qerror(q0.conjugate()*q1);
	vec3_t axis   = qerror.axis ();
	real_t theta  = qerror.angle();
	while(theta >  _pi)
		theta -= 2*_pi;
    while(theta < -_pi)
        theta += 2*_pi;
	
    vec3_t w = theta*axis;
	return (1.0/2.0)*(q0*w + q1*w);
}

inline mat3_t rot_jacobian(vec3_t omega){
	real_t theta = omega.norm();
	if(theta < eps)
		return mat3_t::Identity();

	vec3_t eta = omega/theta;

	mat3_t A = (sin(theta)/theta)*mat3_t::Identity() + (1.0 - (sin(theta)/theta))*vvtrmat(eta, eta) + ((cos(theta) - 1.0)/theta)*cross_mat(eta);
	
	return A;
}

template<class V>
inline V interpolate_pos_linear_diff(real_t t, real_t t0, V p0, real_t t1, V p1){
	real_t h = t1 - t0;
	if(h < eps)
		return p0;
		
	real_t s = (t - t0)/h;
	return (1.0 - s)*p0 + s*p1;
}

template<class V>
inline V interpolate_vel_linear_diff(real_t t0, V p0, real_t t1, V p1){
	real_t h = t1 - t0;
	if(h < eps)
		return V();

	return (p1 - p0)/h;
}

inline quat_t interpolate_slerp_int(real_t t, real_t t0, quat_t q0, vec3_t w0){
	real_t w0norm = w0.norm();
	if(w0norm == 0.0)
		return q0;
	vec3_t axis = w0/w0norm;

	return AngleAxisd(w0norm*(t - t0), axis)*q0;
}

inline quat_t interpolate_slerp_diff(real_t t, real_t t0, quat_t q0, real_t t1, quat_t q1){
    real_t h = t1 - t0;
	if(h < eps)
		return q0;

    real_t s = (t - t0)/h;
    
    AngleAxisd qrel(q0.conjugate()*q1);

	return q0*AngleAxisd(s*qrel.angle(), qrel.axis());
}

inline vec3_t interpolate_angvel_diff(real_t t0, quat_t q0, real_t t1, quat_t q1){
	real_t h = t1 - t0;
	if(h < eps)
		return vec3_t::Zero();

	AngleAxisd qrel(q0.conjugate()*q1);
	vec3_t w = (qrel.angle()/h)*qrel.axis();
	return q0*w;
}

template<class P, class V, class T>
class curve_t{
public:
	typedef P	pos_t;
	typedef V	vel_t;
	typedef T	real_t;

public:
	struct point_t{
		real_t	t;
		pos_t	pos;
		vel_t	vel;

		point_t(real_t _t, pos_t _p, vel_t _v):t(_t), pos(_p), vel(_v){}
	};

	int	type;

	std::vector<struct point_t>	points;

public:
	std::pair<int, int>	GetSegment(real_t t)const{
		if(points.size() < 2)
			return std::make_pair(0, 0);
		int idx = 0;
		while(idx < (int)points.size()-1 && points[idx+1].t < t)
			idx++;
		if(idx == points.size()-1)
			return std::make_pair(idx, idx);
		return std::make_pair(idx, idx+1);
	}
};

template<class V, class T>
class curve_euclid_t : public curve_t<V, V, T>{
public:
    typedef curve_t<V, V, T>  base_t;
	typedef V	pos_t;
	typedef V	vel_t;
	typedef T	real_t;

	pos_t	CalcPos(real_t t){
		std::pair<int,int> seg = this->GetSegment(t);
		struct base_t::point_t& p0 = this->points[seg.first ];
		struct base_t::point_t& p1 = this->points[seg.second];
		return interpolate_pos_linear_diff(t, p0.t, p0.pos, p1.t, p1.pos);
	}

	vel_t	CalcVel(real_t t){
		std::pair<int,int> seg = this->GetSegment(t);
		struct base_t::point_t& p0 = this->points[seg.first ];
		struct base_t::point_t& p1 = this->points[seg.second];
		return interpolate_vel_linear_diff(p0.t, p0.pos, p1.t, p1.pos);
	}

	curve_euclid_t(){
	}
};

class curve_quat_t : public curve_t<quat_t, vec3_t, real_t>{
public:
    typedef typename curve_t<quat_t, vec3_t, real_t> base_t;
	typedef typename curve_t<quat_t, vec3_t, real_t>::point_t point_t;

	quat_t CalcPos(real_t t){
		std::pair<int,int> seg = GetSegment(t);
		point_t& p0 = points[seg.first ];
		point_t& p1 = points[seg.second];
		return interpolate_slerp_diff(t, p0.t, p0.pos, p1.t, p1.pos);
	}

	vel_t	CalcVel(real_t t){
		std::pair<int,int> seg = GetSegment(t);
		point_t& p0 = points[seg.first ];
		point_t& p1 = points[seg.second];
		return interpolate_angvel_diff(p0.t, p0.pos, p1.t, p1.pos);
	}

	curve_quat_t(){
	}
};

typedef curve_euclid_t<real_t, real_t>  curve_real_t;
typedef curve_euclid_t<vec2_t, real_t>  curve_vec2_t;
typedef curve_euclid_t<vec3_t, real_t>  curve_vec3_t;

}
