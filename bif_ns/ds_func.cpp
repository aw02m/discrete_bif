#include "ds_func.hpp"
#include "ds_derivatives.hpp"
#include "dynamical_system.hpp"
#include "eigensolver.hpp"

void store_state(const Eigen::VectorXd &vp, dynamical_system &ds) {
  Eigen::VectorXd x = vp(Eigen::seqN(0, ds.xdim)).real();
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);

  ds.xk[0] = x;
  for (int i = 0; i < ds.period; i++) {
    ds.xk[i + 1] = T(ds.xk[i], ds);
    ds.dTdx[i] = dTdx(ds.xk[i], ds);
    ds.dTdlambda[i] = dTdlambda(ds.xk[i], ds);
    ds.dTdxdx[i] = dTdxdx(ds.xk[i], ds);
    ds.dTdxdlambda[i] = dTdxdlambda(ds.xk[i], ds);
  }

  ds.dTldx = dTldx(ds);
  ds.dTldlambda = dTldlambda(ds);
  if (ds.numerical_diff) {
    ds.dTldxdx = dTldxdx_numeric(ds);
    ds.dTldxdlambda = dTldxdlambda_numeric(ds);
  } else {
    ds.dTldxdx = dTldxdx(ds);
    ds.dTldxdlambda = dTldxdlambda(ds);
  }

  ds.eigvals = eigenvalues(ds);
  ds.mu = Eigen::dcomplex(std::cos(ds.theta), std::sin(ds.theta));

  ds.chara_poly = ds.dTldx - ds.mu * I;
}

Eigen::VectorXd func_newton(const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim + 2);
  Eigen::dcomplex chi(0, 0);

  ret(Eigen::seqN(0, ds.xdim)) = ds.xk[ds.period] - ds.xk[0];
  chi = (ds.chara_poly).determinant();
  ret(ds.xdim) = chi.real();
  ret(ds.xdim + 1) = chi.imag();

  return ret;
}

Eigen::MatrixXd jac_newton(const dynamical_system &ds) {
  Eigen::MatrixXd ret(ds.xdim + 2, ds.xdim + 2);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
  Eigen::MatrixXcd dchidx(1, ds.xdim);
  Eigen::dcomplex dchidlambda(0, 0);
  Eigen::dcomplex dchidtheta(0, 0);
  for (int i = 0; i < ds.xdim; i++) {
    dchidx(0, i) = det_derivative(ds.chara_poly, ds.dTldxdx[i], ds);
  }
  dchidlambda = det_derivative(ds.chara_poly, ds.dTldxdlambda, ds);
  Eigen::MatrixXcd dpolydtheta =
      Eigen::dcomplex(std::sin(ds.theta), -std::cos(ds.theta)) * I;
  dchidtheta = det_derivative(ds.chara_poly, dpolydtheta, ds);

  ret(Eigen::seqN(0, ds.xdim), Eigen::seqN(0, ds.xdim)) = ds.dTldx - I;
  ret(Eigen::seqN(0, ds.xdim), ds.xdim) = ds.dTldlambda;
  ret(Eigen::seqN(0, ds.xdim), ds.xdim + 1) = Eigen::VectorXd::Zero(ds.xdim);
  ret(ds.xdim, Eigen::seqN(0, ds.xdim)) = dchidx.real();
  ret(ds.xdim, ds.xdim) = dchidlambda.real();
  ret(ds.xdim, ds.xdim + 1) = dchidtheta.real();
  ret(ds.xdim + 1, Eigen::seqN(0, ds.xdim)) = dchidx.imag();
  ret(ds.xdim + 1, ds.xdim) = dchidlambda.imag();
  ret(ds.xdim + 1, ds.xdim + 1) = dchidtheta.imag();

  return ret;
}

Eigen::MatrixXd dTldx(const dynamical_system &ds) {
  Eigen::MatrixXd ret = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);

  for (int i = ds.period - 1; i >= 0; i--) {
    ret *= ds.dTdx[i];
  }

  return ret;
}

Eigen::VectorXd dTldlambda(dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim);

  // ds.dTkdlambda is used for dTldxdlambda.
  ret = ds.dTdlambda[0];
  ds.dTkdlambda[0] = ret;
  for (int i = 1; i < ds.period; i++) {
    ret = ds.dTdx[i] * ret + ds.dTdlambda[i];
    ds.dTkdlambda[i] = ret;
  }

  return ret;
}

std::vector<Eigen::MatrixXd> dTldxdx(dynamical_system &ds) {
  std::vector<Eigen::MatrixXd> ret = std::vector<Eigen::MatrixXd>(
      ds.xdim, Eigen::MatrixXd::Zero(ds.xdim, ds.xdim));

  // prods are also used for dTldxdlambda.
  // initialize with Identity matrix.
  ds.frwd_prod = std::vector<Eigen::MatrixXd>(
      ds.period, Eigen::MatrixXd::Identity(ds.xdim, ds.xdim));
  ds.bkwd_prod = std::vector<Eigen::MatrixXd>(
      ds.period, Eigen::MatrixXd::Identity(ds.xdim, ds.xdim));

  for (int i = 0; i < ds.period; i++) {
    for (int j = ds.period - 1; j >= 0; j--) {
      if (j > i) {
        ds.frwd_prod[i] *= ds.dTdx[j];
      } else if (j < i) {
        ds.bkwd_prod[i] *= ds.dTdx[j];
      }
    }
  }

  for (int i = 0; i < ds.xdim; i++) {
    for (int j = 0; j < ds.period; j++) {
      ret[i] +=
          ds.frwd_prod[j] * ds.dTdxdx[j][i] * ds.bkwd_prod[j] * ds.bkwd_prod[j];
    }
  }

  return ret;
}

Eigen::MatrixXd dTldxdlambda(const dynamical_system &ds) {
  Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);
  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);

  ret = ds.dTdxdlambda[0];
  for (int i = 1; i < ds.period; i++) {
    for (int j = 0; j < ds.xdim; j++) {
      temp.col(j) = ds.dTdxdx[i][j] * ds.bkwd_prod[i] * ds.dTkdlambda[i];
    }
    ret = temp + ds.dTdx[i] * ret + ds.dTdxdlambda[i] * ds.bkwd_prod[i];
  }

  return ret;
}

Eigen::dcomplex det_derivative(const Eigen::MatrixXcd &A,
                               const Eigen::MatrixXcd &dA,
                               const dynamical_system &ds) {
  Eigen::MatrixXcd temp(ds.xdim, ds.xdim);
  Eigen::dcomplex ret(0, 0);

  for (int i = 0; i < ds.xdim; i++) {
    temp = A;
    temp.col(i) = dA.col(i).cast<Eigen::dcomplex>();
    ret += temp.determinant();
  }

  return ret;
}

// numerical differentiation version
std::vector<Eigen::MatrixXd> dTldxdx_numeric(dynamical_system &ds) {
  std::vector<Eigen::MatrixXd> ret = std::vector<Eigen::MatrixXd>(
      ds.xdim, Eigen::MatrixXd::Zero(ds.xdim, ds.xdim));
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);
  Eigen::MatrixXd Ah = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
  Eigen::VectorXd h_vec(ds.xdim);
  double h = ds.dif_strip;

  for (int i = 0; i < ds.xdim; i++) {
    for (int j = 0; j < ds.xdim; j++) {
      if (i == j) {
        h_vec(j) = h;
      } else {
        h_vec(j) = 0;
      }
    }
    A = ds.dTldx;
    Ah = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
    for (int j = ds.period - 1; j >= 0; j--) {
      Ah *= dTdx(ds.xk[j] + h_vec, ds);
    }
    ret[i] = (Ah - A) / h;
  }

  return ret;
}

Eigen::MatrixXd dTldxdlambda_numeric(dynamical_system &ds) {
  Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);
  Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);
  Eigen::MatrixXd A = ds.dTldx;
  Eigen::MatrixXd Ah = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
  double lambda = ds.params(ds.var_param);
  double h = ds.dif_strip;

  ds.params(ds.var_param) += h;
  for (int i = ds.period - 1; i >= 0; i--) {
    Ah *= dTdx(ds.xk[i], ds);
  }
  ret = (Ah - A) / h;
  ds.params(ds.var_param) = lambda;

  return ret;
}