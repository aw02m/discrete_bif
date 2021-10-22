/* GENERALIZED HENON MAP */

#include "dynamical_system.hpp"

void dynamical_system::function(const Eigen::VectorXd &x) {
  /* Zero Initialization */
  T = Eigen::VectorXd::Zero(xdim);
  dTdx = Eigen::MatrixXd::Zero(xdim, xdim);
  dTdlambda = Eigen::VectorXd::Zero(xdim);
  dTdxdlambda = Eigen::MatrixXd::Zero(xdim, xdim);
  dTdxdx =
      std::vector<Eigen::MatrixXd>(xdim, Eigen::MatrixXd::Zero(xdim, xdim));
  /* End Zero Initialization */

  static double a, b, c, d;
  static unsigned int m;
  a = p(0);
  b = p(1);
  c = p(2);
  d = p(3);
  m = p(4);

  /* P(x) version */
  // T(0) = d * x(1);
  // T(1) = a * std::pow(x(0), m) * (x(0) * x(0) - b * b) + c * x(0);

  // dTdx(0, 1) = d;
  // dTdx(1, 0) = a * m * std::pow(x(0), m - 1) * (x(0) * x(0) - b * b) +
  //              a * std::pow(x(0), m) * 2 * x(0) + c;

  // switch (var_param) {
  // case 0:
  //   dTdlambda(1) = std::pow(x(0), m) * (x(0) * x(0) - b * b);
  //   break;
  // case 3:
  //   dTdlambda(0) = x(1);
  //   break;
  // }

  // dTdxdx[0](1, 0) =
  //     a * m * (m - 1) * std::pow(x(0), m - 2) * (x(0) * x(0) - b * b) +
  //     a * m * std::pow(x(0), m - 1) * 2 * x(0) +
  //     a * m * std::pow(x(0), m - 1) * 2 * x(0) + a * std::pow(x(0), m) * 2;

  // switch (var_param) {
  // case 0:
  //   dTdxdlambda(1, 0) = m * std::pow(x(0), m - 1) * (x(0) * x(0) - b * b) +
  //                       2 * std::pow(x(0), m + 1);
  //   break;
  // case 3:
  //   dTdxdlambda(0, 1) = 1;
  //   break;
  // }

  /* P(y) version */
  T(0) = d * x(1);
  T(1) = a * std::pow(x(1), m) * (x(1) * x(1) - b * b) + c * x(0);

  dTdx(0, 1) = d;
  dTdx(1, 0) = c;
  dTdx(1, 1) = a * m * std::pow(x(1), m - 1) * (x(1) * x(1) - b * b) +
               a * std::pow(x(1), m) * 2 * x(1);

  switch (var_param) {
  case 0:
    dTdlambda(1) = std::pow(x(1), m) * (x(1) * x(1) - b * b);
    break;
  case 3:
    dTdlambda(0) = x(1);
    break;
  }

  dTdxdx[1](1, 1) =
      a * m * (m - 1) * std::pow(x(1), m - 2) * (x(1) * x(1) - b * b) +
      2 * a * m * std::pow(x(1), m) + 2 * a * (m + 1) * std::pow(x(1), m);

  switch (var_param) {
  case 0:
    dTdxdlambda(1, 1) = m * std::pow(x(1), m - 1) * (x(1) * x(1) - b * b) +
                        2 * std::pow(x(1), m + 1);
    break;
  case 3:
    dTdxdlambda(0, 1) = 1;
    break;
  }
}