#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

#include "sys_common.hpp"

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
  bool numerical_diff;
  double dif_strip;
  double eps;
  double explode;

  Eigen::VectorXd x0;
  Eigen::VectorXd params;

  Eigen::VectorXcd eigvals;
  Eigen::dcomplex mu;
  double theta;

  std::vector<Eigen::VectorXd> xk;

  std::vector<Eigen::MatrixXd> dTdx;
  std::vector<Eigen::VectorXd> dTdlambda;
  std::vector<std::vector<Eigen::MatrixXd>> dTdxdx;
  std::vector<Eigen::MatrixXd> dTdxdlambda;

  // these are common terms of second derivatives.
  std::vector<Eigen::MatrixXd> frwd_prod;
  std::vector<Eigen::MatrixXd> bkwd_prod;
  std::vector<Eigen::VectorXd> dTkdlambda;

  Eigen::MatrixXd dTldx;
  Eigen::MatrixXd dTldlambda;
  std::vector<Eigen::MatrixXd> dTldxdx;
  Eigen::MatrixXd dTldxdlambda;

  Eigen::MatrixXcd chara_poly;
};

#endif