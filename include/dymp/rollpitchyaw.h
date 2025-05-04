#pragma once

#include <dymp/types.h>

namespace dymp{

    /** 
     *  @file
     *  conversion between roll-pitch-yaw and quaternion
     *  q = Rz(yaw)*Ry(pitch)*Rx(roll)
     **/

    /// quaternion to roll-pitch-yaw
    vec3_t ToRollPitchYaw(const quat_t& q);

    /// roll-pitch-yaw to quaternion
    quat_t FromRollPitchYaw(const vec3_t& angles);

    // time derivative of roll-pitch-yaw to angular velocity
    vec3_t VelocityFromRollPitchYaw(const vec3_t& angle, const vec3_t& angled);

    // jacobian from rpy derivative to angular velocity
    mat3_t JacobianFromRollPitchYaw(const vec3_t& angle);

}