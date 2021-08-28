#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

#include "sys_common.hpp"
#include <nlohmann/json.hpp>

class dynamical_system {
public:
  dynamical_system(nlohmann::json json);

  unsigned int xdim;
  unsigned int period;

  // Newton param
  unsigned int inc_param;
  unsigned int var_param;
  double delta_inc;
  unsigned int inc_iter;
  unsigned int max_iter;
  double dif_strip;
  double eps;
  double explode;

  Eigen::VectorXd x0;
  Eigen::VectorXd params;

  std::vector<Eigen::VectorXd> xk;
  std::vector<Eigen::VectorXd> fk;

  std::vector<Eigen::MatrixXd> dfdx;

  std::vector<Eigen::MatrixXd> dTdx;
  std::vector<std::vector<Eigen::MatrixXd>> dTdxdx;
  std::vector<Eigen::VectorXd> dTdlambda;
  std::vector<Eigen::MatrixXd> dTdxdlambda;

  Eigen::MatrixXd dTldx;

  Eigen::VectorXcd eigvals;
};

#endif