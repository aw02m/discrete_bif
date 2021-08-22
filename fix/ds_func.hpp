#ifndef DS_FUNC_HPP_
#define DS_FUNC_HPP_

#include "sys_common.hpp"

class dynamical_system;

void store_state(const Eigen::VectorXd &vp, dynamical_system &ds);

Eigen::VectorXd func_newton(const dynamical_system &ds);
Eigen::MatrixXd jac_newton(const dynamical_system &ds);

Eigen::MatrixXd dTldx(const dynamical_system &ds);

#endif