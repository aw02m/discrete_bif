#include "dynamical_system.hpp"

dynamical_system::dynamical_system(nlohmann::json json) {
  xdim = json["x0"].size();
  period = json["period"];
  inc_param = json["inc_param"];
  var_param = json["var_param"];
  delta_inc = json["delta_inc"];
  inc_iter = json["inc_iter"];
  max_iter = json["max_iter"];
  dif_strip = json["dif_strip"];
  eps = json["eps"];
  explode = json["explode"];

  /* These json array should be casted to the STL container type*/
  std::vector<double> fixed_arr = json["x0"];
  Eigen::Map<Eigen::VectorXd> x0(fixed_arr.data(), fixed_arr.size());
  this->x0 = x0;

  std::vector<double> params_arr = json["params"];
  Eigen::Map<Eigen::VectorXd> params(params_arr.data(), params_arr.size());
  this->params = params;

  xk = std::vector<Eigen::VectorXd>(period + 1, Eigen::VectorXd::Zero(xdim));
  dTdx =
      std::vector<Eigen::MatrixXd>(period, Eigen::MatrixXd::Zero(xdim, xdim));
  dTdlambda = std::vector<Eigen::VectorXd>(period, Eigen::VectorXd::Zero(xdim));
  dTdxdx = std::vector<std::vector<Eigen::MatrixXd>>(
      period,
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim)));
  dTdxdlambda =
      std::vector<Eigen::MatrixXd>(period, Eigen::MatrixXd::Zero(xdim, xdim));
  dTldxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));

  frwd_prod = std::vector<Eigen::MatrixXd>(
      period, Eigen::MatrixXd::Identity(xdim, xdim));
  bkwd_prod = std::vector<Eigen::MatrixXd>(
      period, Eigen::MatrixXd::Identity(xdim, xdim));
  dTkdlambda =
      std::vector<Eigen::VectorXd>(period, Eigen::VectorXd::Zero(xdim));
}