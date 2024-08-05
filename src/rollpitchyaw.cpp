#include <rollpitchyaw.h>
#include <util.h>

using namespace Eigen;

namespace dymp{

vec3_t ToRollPitchYaw(const quat_t& q){
	vec3_t angles;

	vec3_t xdir = q*vec3_t::UnitX();
	angles[2] = atan2( xdir.y(), xdir.x());
	angles[1] = atan2(-xdir.z(), sqrt(xdir.x()*xdir.x() + xdir.y()*xdir.y()));
	
	quat_t qroll = AngleAxisd(-angles[1], vec3_t::UnitY()) * AngleAxisd(-angles[2], vec3_t::UnitZ()) * q;
	angles[0] = 2.0 * atan2(qroll.x(), qroll.w());

	// yaw angle needs wrapping
	if(angles[0] >  _pi) angles[0] -= 2.0*_pi;
	if(angles[0] < -_pi) angles[0] += 2.0*_pi;

	return angles;
}

quat_t FromRollPitchYaw(const vec3_t& angles){
	return AngleAxisd(angles[2], vec3_t::UnitZ())
         * AngleAxisd(angles[1], vec3_t::UnitY())
         * AngleAxisd(angles[0], vec3_t::UnitX());
}

vec3_t VelocityFromRollPitchYaw(const vec3_t& angle, const vec3_t& angled){
    AngleAxisd qz  = AngleAxisd(angle.z(), ez);
	quat_t     qzy = qz*AngleAxisd(angle.y(), ey);
	vec3_t wx = ex*angled.x();
	vec3_t wy = ey*angled.y();
	vec3_t wz = ez*angled.z();

	return wz + qz*wy + qzy*wx;
}

}