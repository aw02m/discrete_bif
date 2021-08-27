#ifndef DS_FUNC_HPP_
#define DS_FUNC_HPP_

#include "sys_common.hpp"

class dynamical_system;

void store_state(const Eigen::VectorXd &vp, dynamical_system &ds);

Eigen::VectorXd func_newton(const dynamical_system &ds);
Eigen::MatrixXd jac_newton(const dynamical_system &ds);

Eigen::MatrixXd dTldx(const dynamical_system &ds);
Eigen::VectorXd dTldlambda(dynamical_system &ds);
std::vector<Eigen::MatrixXd> dTldxdx(dynamical_system &ds);
Eigen::MatrixXd dTldxdlambda(const dynamical_system &ds);
// numerical differentiation version
std::vector<Eigen::MatrixXd> dTldxdx_numeric(dynamical_system &ds);
Eigen::MatrixXd dTldxdlambda_numeric(dynamical_system &ds);

Eigen::dcomplex det_derivative(const Eigen::MatrixXcd &A,
                               const Eigen::MatrixXcd &dA,
                               const dynamical_system &ds);

#endif