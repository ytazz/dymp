#pragma once

#include <memory>

#include <Eigen/Eigen>

namespace dymp{;

typedef double                     real_t;
typedef Eigen::Vector2d            vec2_t;
typedef Eigen::Vector3d            vec3_t;
typedef Eigen::Quaterniond         quat_t;
typedef Eigen::Matrix2d            mat2_t;
typedef Eigen::Matrix3d            mat3_t;
typedef Eigen::Matrix<double,6,6>  mat6_t;
typedef Eigen::Matrix<double,6,3>  mat63_t;
typedef Eigen::Matrix<double,3,6>  mat36_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1>  vvec_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  vmat_t;

}
