#include "ds_func.hpp"

void store_state(const Eigen::VectorXd &vp, dynamical_system &ds) {
  Eigen::VectorXd x = vp(Eigen::seqN(0, ds.xdim));
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);

  ds.xk[0] = x;
  for (int i = 0; i < ds.period; i++) {
    ds.xk[i + 1] = T(ds.xk[i], ds);
    ds.dTdx[i] = dTdx(ds.xk[i], ds);
  }

  ds.dTldx = dTldx(ds);
}

Eigen::VectorXd func_newton(const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim);

  ret(Eigen::seqN(0, ds.xdim)) = ds.xk[ds.period] - ds.xk[0];

  return ret;
}

Eigen::MatrixXd jac_newton(const dynamical_system &ds) {
  Eigen::MatrixXd ret(ds.xdim, ds.xdim);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);

  ret(Eigen::seqN(0, ds.xdim), Eigen::seqN(0, ds.xdim)) = ds.dTldx - I;

  return ret;
}

Eigen::MatrixXd dTldx(const dynamical_system &ds) {
  Eigen::MatrixXd ret = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);

  for (int i = ds.period - 1; i >= 0; i--) {
    ret *= ds.dTdx[i];
  }

  return ret;
}