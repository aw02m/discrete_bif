#ifndef DS_DERIVATIVES_HPP_
#define DS_DERIVATIVES_HPP_

#include "sys_common.hpp"

class dynamical_system;

Eigen::VectorXd T(const Eigen::VectorXd &x, const dynamical_system &ds);
Eigen::MatrixXd dTdx(const Eigen::VectorXd &x, const dynamical_system &ds);
Eigen::VectorXd dTdlambda(const Eigen::VectorXd &x, const dynamical_system &ds);
std::vector<Eigen::MatrixXd> dTdxdx(const Eigen::VectorXd &x,
                                    const dynamical_system &ds);
Eigen::MatrixXd dTdxdlambda(const Eigen::VectorXd &x,
                            const dynamical_system &ds);

#endif