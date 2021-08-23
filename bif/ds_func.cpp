#include "ds_func.hpp"

void store_state(const Eigen::VectorXd &vp, dynamical_system &ds) {
  Eigen::VectorXd x = vp(Eigen::seqN(0, ds.xdim));
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
  ds.dTldxdx = dTldxdx(ds);
  ds.dTldxdlambda = dTldxdlambda(ds);
}

Eigen::VectorXd func_newton(const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim + 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
  double mu = -1;

  ret(Eigen::seqN(0, ds.xdim)) = ds.xk[ds.period] - ds.xk[0];
  ret(ds.xdim) = (ds.dTldx - mu * I).determinant();

  return ret;
}

Eigen::MatrixXd jac_newton(const dynamical_system &ds) {
  Eigen::MatrixXd ret(ds.xdim + 1, ds.xdim + 1);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
  double mu = -1;
  Eigen::MatrixXd A = ds.dTldx - mu * I;

  ret(Eigen::seqN(0, ds.xdim), Eigen::seqN(0, ds.xdim)) = ds.dTldx - I;
  ret(Eigen::seqN(0, ds.xdim), ds.xdim) = ds.dTldlambda;
  for (int i = 0; i < ds.xdim; i++) {
    ret(ds.xdim, i) = det_derivative(A, ds.dTldxdx[i], ds);
  }
  ret(ds.xdim, ds.xdim) = det_derivative(A, ds.dTldxdlambda, ds);

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

  // temp.col(0) = ds.dTdxdx[1][0] * ds.dTdlambda[0];
  // temp.col(1) = ds.dTdxdx[1][1] * ds.dTdlambda[0];

  // ret = temp + ds.dTdxdlambda[1];
  // ret *= ds.dTdx[0];
  // ret += ds.dTdx[1] * ds.dTdxdlambda[0];

  return ret;
}

double det_derivative(const Eigen::MatrixXd &A, const Eigen::MatrixXd &dA,
                      const dynamical_system &ds) {
  Eigen::MatrixXd temp(ds.xdim, ds.xdim);
  double ret = 0;

  for (int i = 0; i < ds.xdim; i++) {
    temp = A;
    temp.col(i) = dA.col(i);
    ret += temp.determinant();
  }

  return ret;
}

// numerical differentiation version
// std::vector<Eigen::MatrixXd> dTldxdx(dynamical_system &ds) {
//   std::vector<Eigen::MatrixXd> ret = std::vector<Eigen::MatrixXd>(
//       ds.xdim, Eigen::MatrixXd::Zero(ds.xdim, ds.xdim));
//   Eigen::MatrixXd A = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);
//   Eigen::MatrixXd Ah = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
//   Eigen::VectorXd h_vec(ds.xdim);
//   double h = ds.dif_strip;

//   for (int i = 0; i < ds.xdim; i++) {
//     for (int j = 0; j < ds.xdim; j++) {
//       if (i == j) {
//         h_vec(j) = h;
//       } else {
//         h_vec(j) = 0;
//       }
//     }
//     A = ds.dTldx;
//     Ah = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
//     for (int j = ds.period - 1; j >= 0; j--) {
//       Ah *= dTdx(ds.xk[j] + h_vec, ds);
//     }
//     ret[i] = (Ah - A) / h;
//   }

//   return ret;
// }

// Eigen::MatrixXd dTldxdlambda(dynamical_system &ds) {
//   Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);
//   Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(ds.xdim, ds.xdim);
//   Eigen::MatrixXd A = ds.dTldx;
//   Eigen::MatrixXd Ah = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
//   double lambda = ds.params(ds.var_param);
//   double h = ds.dif_strip;

//   ds.params(ds.var_param) += h;
//   for (int i = ds.period - 1; i >= 0; i--) {
//     Ah *= dTdx(ds.xk[i], ds);
//   }
//   ret = (Ah - A) / h;
//   ds.params(ds.var_param) = lambda;

//   return ret;
// }