#include "dynamical_system.hpp"

Eigen::VectorXd dynamical_system::function(const Eigen::VectorXd &x) {
  Eigen::VectorXd f(xdim);

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

  /* P(y) version */
  f(0) = d * x(1);
  f(1) = a * std::pow(x(1), m) * (x(1) * x(1) - b * b) + c * x(0);

  return f;
}