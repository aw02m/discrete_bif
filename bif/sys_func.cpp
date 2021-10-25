/* GENERALIZED HENON MAP */

#include "dynamical_system.hpp"

void dynamical_system::function([[maybe_unused]] int k, Eigen::VectorXd &x) {
  /* Zero Initialization */
  f = Eigen::VectorXd::Zero(xdim);
  dfdx = Eigen::MatrixXd::Zero(xdim, xdim);
  dfdlambda = Eigen::VectorXd::Zero(xdim);
  dfdxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));
  dfdxdlambda = Eigen::MatrixXd::Zero(xdim, xdim);
  /* End Zero Initialization */

  // parameters
  static double a, b, c, d;
  static unsigned int m;
  a = p(0);
  b = p(1);
  c = p(2);
  d = p(3);
  m = p(4);

  /* P(y) version */
  f(0) = d * x(1);
  f(1) = a * std::pow(x(1), m) * (x(1) * x(1) - b * b) + c * x(0);

  dfdx(0, 1) = d;
  dfdx(1, 0) = c;
  dfdx(1, 1) = a * m * std::pow(x(1), m - 1) * (x(1) * x(1) - b * b) +
               a * std::pow(x(1), m) * 2 * x(1);

  switch (var_param) {
  case 0:
    dfdlambda(1) = std::pow(x(1), m) * (x(1) * x(1) - b * b);
    break;
  case 3:
    dfdlambda(0) = x(1);
    break;
  }

  // numerical_diff does not require 2nd derivatives.
  if (!numerical_diff) {
    dfdxdx[1](1, 1) =
        a * m * (m - 1) * std::pow(x(1), m - 2) * (x(1) * x(1) - b * b) +
        2 * a * m * std::pow(x(1), m) + 2 * a * (m + 1) * std::pow(x(1), m);

    switch (var_param) {
    case 0:
      dfdxdlambda(1, 1) = m * std::pow(x(1), m - 1) * (x(1) * x(1) - b * b) +
                          2 * std::pow(x(1), m + 1);
      break;
    case 3:
      dfdxdlambda(0, 1) = 1;
      break;
    }
  }

  /******************************************************/
  /*     DO NOT EDIT BELOW (variational equations)      */
  /******************************************************/
  unsigned int counter = xdim;
  // variational state (transform to matrix shape for easy producting)
  Eigen::MatrixXd state_dphidx = x(Eigen::seqN(counter, size_dphidx));
  state_dphidx.resize(xdim, xdim);
  counter += size_dphidx;
  Eigen::VectorXd state_dphidlambda = x(Eigen::seqN(counter, size_dphidlambda));
  counter += size_dphidlambda;
  std::vector<Eigen::MatrixXd> state_dphidxdx(
      xdim, Eigen::MatrixXd::Zero(xdim, xdim));
  Eigen::MatrixXd temp;
  for (int i = 0; i < xdim; i++) {
    temp = x(Eigen::seqN(counter + size_dphidx * i, size_dphidx));
    temp.resize(xdim, xdim);
    state_dphidxdx[i] = temp;
    temp.resize(size_dphidx, 1);
  }
  counter += size_dphidxdx;
  Eigen::MatrixXd state_dphidxdlambda =
      x(Eigen::seqN(counter, size_dphidxdlambda));

  counter = 0;

  // phi
  x(Eigen::seqN(counter, xdim)) = f;
  counter += xdim;

  // dphidx
  x(Eigen::seqN(counter, size_dphidx)) = dfdx * state_dphidx;
  counter += size_dphidx;

  // dphidlambda
  x(Eigen::seqN(counter, xdim)) = dfdx * state_dphidlambda + dfdlambda;
  counter += size_dphidlambda;

  // numerical_diff does not require 2nd derivatives.
  if (!numerical_diff) {
    // dphidxdx
    for (int i = 0; i < xdim; i++) {
      x(Eigen::seqN(counter, size_dphidx)) =
          dfdxdx[i] * state_dphidx * state_dphidx + dfdx * state_dphidxdx[i];
      counter += size_dphidx;
    }

    // dphidxdlambda
    for (int i = 0; i < xdim; i++) {
      x(Eigen::seqN(counter, xdim)) =
          dfdxdx[i] * state_dphidlambda * state_dphidx +
          dfdxdlambda * state_dphidx.col(i) + dfdx * state_dphidxdlambda.col(i);
      counter += xdim;
    }
  }
}