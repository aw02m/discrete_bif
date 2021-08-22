#include "ds_derivatives.hpp"

Eigen::VectorXd T(const Eigen::VectorXd &x, const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim);
  double a, b, c, d;
  unsigned int m;
  a = ds.params(0);
  b = ds.params(1);
  c = ds.params(2);
  d = ds.params(3);
  m = ds.params(4);

  ret(0) = d * x(1);
  ret(1) = a * std::pow(x(1), m) * (x(1) * x(1) - b * b) + c * x(0);

  return ret;
}

Eigen::MatrixXd dTdx(const Eigen::VectorXd &x, const dynamical_system &ds) {
  Eigen::MatrixXd ret(ds.xdim, ds.xdim);
  double a, b, c, d;
  unsigned int m;
  a = ds.params(0);
  b = ds.params(1);
  c = ds.params(2);
  d = ds.params(3);
  m = ds.params(4);

  ret(0, 0) = 0;
  ret(0, 1) = d;
  ret(1, 0) = c;
  ret(1, 1) = a * m * std::pow(x(1), m - 1) * (x(1) * x(1) - b * b) +
              a * std::pow(x(1), m) * 2 * x(1);

  return ret;
}

Eigen::VectorXd dTdlambda(const Eigen::VectorXd &x,
                          const dynamical_system &ds) {
  Eigen::VectorXd ret(ds.xdim);
  double b = ds.params(1);
  unsigned int m = ds.params(4);

  switch (ds.var_param) {
  case 0:
    ret(0) = 0;
    ret(1) = std::pow(x(1), m) * (x(1) * x(1) - b * b);
    break;
  case 3:
    ret(0) = x(1);
    ret(1) = 0;
    break;
  }

  return ret;
}

std::vector<Eigen::MatrixXd> dTdxdx(const Eigen::VectorXd &x,
                                    const dynamical_system &ds) {
  std::vector<Eigen::MatrixXd> ret(ds.xdim,
                                   Eigen::MatrixXd::Zero(ds.xdim, ds.xdim));
  double a, b, c, d;
  unsigned int m;
  a = ds.params(0);
  b = ds.params(1);
  c = ds.params(2);
  d = ds.params(3);
  m = ds.params(4);

  ret[0](0, 0) = 0;
  ret[0](0, 1) = 0;
  ret[0](1, 0) = 0;
  ret[0](1, 1) = 0;
  ret[1](0, 0) = 0;
  ret[1](0, 1) = 0;
  ret[1](1, 0) = 0;
  ret[1](1, 1) = a * m * (m - 1) * std::pow(x(1), m - 2) * (x(1) * x(1) - b * b) +
                 2 * a * m * std::pow(x(1), m) +
                 2 * a * (m + 1) * std::pow(x(1), m);

  return ret;
}

Eigen::MatrixXd dTdxdlambda(const Eigen::VectorXd &x,
                            const dynamical_system &ds) {
  Eigen::MatrixXd ret(ds.xdim, ds.xdim);
  double b = ds.params(1);
  unsigned int m = ds.params(4);

  switch (ds.var_param) {
  case 0:
    ret(0, 0) = 0;
    ret(0, 1) = 0;
    ret(1, 0) = 0;
    ret(1, 1) = m * std::pow(x(0), m - 1) * (x(0) * x(0) - b * b) +
                2 * std::pow(x(0), m + 1);
    break;
  case 3:
    ret(0, 0) = 0;
    ret(0, 1) = 1;
    ret(1, 0) = 0;
    ret(1, 1) = 0;
    break;
  }

  return ret;
}

// Eigen::VectorXd T(const Eigen::VectorXd &x, const dynamical_system &ds) {
//   Eigen::VectorXd ret(ds.xdim);
//   double a, b, c, d;
//   unsigned int m;
//   a = ds.params(0);
//   b = ds.params(1);
//   c = ds.params(2);
//   d = ds.params(3);
//   m = ds.params(4);

//   ret(0) = d * x(1);
//   ret(1) = a * std::pow(x(0), m) * (x(0) * x(0) - b * b) + c * x(0);

//   return ret;
// }

// Eigen::MatrixXd dTdx(const Eigen::VectorXd &x, const dynamical_system &ds) {
//   Eigen::MatrixXd ret(ds.xdim, ds.xdim);
//   double a, b, c, d;
//   unsigned int m;
//   a = ds.params(0);
//   b = ds.params(1);
//   c = ds.params(2);
//   d = ds.params(3);
//   m = ds.params(4);

//   ret(0, 0) = 0;
//   ret(0, 1) = d;
//   ret(1, 0) = a * m * std::pow(x(0), m - 1) * (x(0) * x(0) - b * b) +
//               a * std::pow(x(0), m) * 2 * x(0) + c;
//   ret(1, 1) = 0;

//   return ret;
// }

// Eigen::VectorXd dTdlambda(const Eigen::VectorXd &x,
//                           const dynamical_system &ds) {
//   Eigen::VectorXd ret(ds.xdim);
//   double b = ds.params(1);
//   unsigned int m = ds.params(4);

//   switch (ds.var_param) {
//   case 0:
//     ret(0) = 0;
//     ret(1) = std::pow(x(0), m) * (x(0) * x(0) - b * b);
//     break;
//   case 3:
//     ret(0) = x(1);
//     ret(1) = 0;
//     break;
//   }

//   return ret;
// }

// std::vector<Eigen::MatrixXd> dTdxdx(const Eigen::VectorXd &x,
//                                     const dynamical_system &ds) {
//   std::vector<Eigen::MatrixXd> ret(ds.xdim,
//                                    Eigen::MatrixXd::Zero(ds.xdim, ds.xdim));
//   double a, b, c, d;
//   unsigned int m;
//   a = ds.params(0);
//   b = ds.params(1);
//   c = ds.params(2);
//   d = ds.params(3);
//   m = ds.params(4);

//   ret[0](1, 0) =
//       a * m * (m - 1) * std::pow(x(0), m - 2) * (x(0) * x(0) - b * b) +
//       a * m * std::pow(x(0), m - 1) * 2 * x(0) +
//       a * m * std::pow(x(0), m - 1) * 2 * x(0) + a * std::pow(x(0), m) * 2;

//   return ret;
// }

// Eigen::MatrixXd dTdxdlambda(const Eigen::VectorXd &x,
//                             const dynamical_system &ds) {
//   Eigen::MatrixXd ret = Eigen::MatrixXd::Identity(ds.xdim, ds.xdim);
//   double b = ds.params(1);
//   unsigned int m = ds.params(4);

//   switch (ds.var_param) {
//   case 0:
//     ret(1, 0) = m * std::pow(x(0), m - 1) * (x(0) * x(0) - b * b) +
//                 2 * std::pow(x(0), m + 1);
//     break;
//   case 3:
//     ret(0, 1) = 1;
//     break;
//   }

//   return ret;
// }