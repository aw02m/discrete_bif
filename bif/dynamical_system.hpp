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
  bool numerical_diff;
  double dif_strip;
  double eps;
  double explode;
  bool fix_mode;

  std::string out_path;

  Eigen::VectorXd x0;
  Eigen::VectorXd params;
  Eigen::VectorXd p;

  Eigen::VectorXcd eigvals;
  Eigen::dcomplex mu;
  double theta;

  Eigen::VectorXd T;
  Eigen::MatrixXd dTdx;
  Eigen::VectorXd dTdlambda;
  std::vector<Eigen::MatrixXd> dTdxdx;
  Eigen::MatrixXd dTdxdlambda;

  std::vector<Eigen::VectorXd> xk;
  std::vector<Eigen::MatrixXd> dTdx_arr;
  std::vector<Eigen::VectorXd> dTdlambda_arr;
  std::vector<std::vector<Eigen::MatrixXd>> dTdxdx_arr;
  std::vector<Eigen::MatrixXd> dTdxdlambda_arr;

  // these are common terms of second derivatives.
  std::vector<Eigen::MatrixXd> frwd_prod_arr;
  std::vector<Eigen::MatrixXd> bkwd_prod_arr;
  std::vector<Eigen::VectorXd> dTkdlambda_arr;

  Eigen::MatrixXd dTldx;
  Eigen::MatrixXd dTldlambda;
  std::vector<Eigen::MatrixXd> dTldxdx;
  Eigen::MatrixXd dTldxdlambda;

  Eigen::MatrixXcd chara_poly;

  void function(const Eigen::VectorXd &x);
  Eigen::dcomplex det_derivative(const Eigen::MatrixXcd &A, const Eigen::MatrixXcd &dA);

  void store_states(const Eigen::VectorXd &v);
  void store_states_numeric(const Eigen::VectorXd &v);
  Eigen::VectorXd newton_F();
  Eigen::MatrixXd newton_J();
  void store_states_fix(const Eigen::VectorXd &v);
  Eigen::VectorXd newton_fix_F();
  Eigen::MatrixXd newton_fix_J();
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
  newton_FJ(const Eigen::VectorXd &v);
};

#endif