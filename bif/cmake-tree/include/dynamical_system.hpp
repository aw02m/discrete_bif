#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

#include "sys_common.hpp"
#include <nlohmann/json.hpp>

class dynamical_system {
public:
  dynamical_system(nlohmann::json json);

  unsigned int xdim;
  unsigned int bialt_dim;
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
  unsigned int mode;

  std::string out_path;
  std::string json_out_path;

  Eigen::VectorXd x0;
  Eigen::VectorXd p;

  Eigen::VectorXcd eigvals;

  Eigen::VectorXd f;
  Eigen::MatrixXd dfdx;
  Eigen::VectorXd dfdlambda;
  std::vector<Eigen::MatrixXd> dfdxdx;
  Eigen::MatrixXd dfdxdlambda;

  std::vector<Eigen::VectorXd> xk;
  unsigned int size_dphidx;
  unsigned int size_dphidlambda;
  unsigned int size_dphidxdx;
  unsigned int size_dphidxdlambda;
  Eigen::MatrixXd dphidx;
  Eigen::VectorXd dphidlambda;
  std::vector<Eigen::MatrixXd> dphidxdx;
  Eigen::MatrixXd dphidxdlambda;

  Eigen::MatrixXd dTldx;
  Eigen::MatrixXd dTldlambda;
  std::vector<Eigen::MatrixXd> dTldxdx;
  Eigen::MatrixXd dTldxdlambda;

  Eigen::MatrixXd chara_poly;

  void function([[maybe_unused]] int k, Eigen::VectorXd &x);

  void store_states(const Eigen::VectorXd &v);
  void store_states_numeric(const Eigen::VectorXd &v);
  Eigen::VectorXd newton_F_G();
  Eigen::MatrixXd newton_J_G();
  Eigen::VectorXd newton_F_PD();
  Eigen::MatrixXd newton_J_PD();
  Eigen::VectorXd newton_F_NS();
  Eigen::MatrixXd newton_J_NS();
  void store_states_fix(const Eigen::VectorXd &v);
  Eigen::VectorXd newton_F_fix();
  Eigen::MatrixXd newton_J_fix();
  std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
  newton_FJ(const Eigen::VectorXd &v);
};

#endif