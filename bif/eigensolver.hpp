#ifndef DS_EIGENSOLVER_HPP_
#define DS_EIGENSOLVER_HPP_

#include "sys_common.hpp"

class dynamical_system;

Eigen::VectorXcd eigenvalues(const dynamical_system &ds);
Eigen::dcomplex bifeigvals(const dynamical_system &ds);

#endif