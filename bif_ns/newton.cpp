#include "newton.hpp"
#include <chrono>

void newton(dynamical_system &ds) {
  Eigen::VectorXd vp(ds.xdim + 2);
  // vp << ds.x0 << ds.params(ds.var_param);
  vp(Eigen::seqN(0, ds.xdim)) = ds.x0;
  vp(ds.xdim) = ds.params(ds.var_param);
  vp(ds.xdim + 1) = ds.theta;
  Eigen::VectorXd vn(ds.xdim + 2);
  Eigen::VectorXd F(ds.xdim + 2);
  Eigen::MatrixXd J(ds.xdim + 2, ds.xdim + 2);
  double norm;

  for (int p = 0; p < ds.inc_iter; p++) {
    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < ds.max_iter; i++) {
      store_state(vp, ds);
      F = func_newton(ds);
      J = jac_newton(ds);
      vn = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(J).solve(-F) + vp;

      norm = (vn - vp).norm();
      if (norm < ds.eps) {
        auto end = std::chrono::system_clock::now();
        auto dur = end - start;
        auto msec =
            std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
        std::cout << "**************************************************"
                  << std::endl;
        std::cout << p << " : converged (iter = " << i + 1 << ", ";
        std::cout << "time = " << msec << "[msec])" << std::endl;
        std::cout << "params : " << ds.params.transpose() << std::endl;
        std::cout << "x0     : "
                  << vn(Eigen::seqN(0, ds.xdim)).transpose()
                  << std::endl;
        std::cout << "(Re(μ), Im(μ)), abs(μ), arg(μ) :" << std::endl;
        for (int k = 0; k < ds.xdim; k++) {
          std::cout << ds.eigvals(k) << ", ";
          std::cout << std::abs(ds.eigvals(k)) << ", ";
          std::cout << std::arg(ds.eigvals(k)) * (180 / EIGEN_PI) << std::endl;
        }
        std::cout << "**************************************************"
                  << std::endl;
        vp = vn;
        ds.params(ds.var_param) = vn(ds.xdim);
        ds.theta = vn(ds.xdim + 1);
        break;
      } else if (norm >= ds.explode) {
        std::cerr << "explode (iter = " << i + 1 << ")" << std::endl;
        exit(1);
      }

      if (i >= ds.max_iter - 1) {
        std::cerr << "iter over (iter = " << i + 1 << ")" << std::endl;
        exit(1);
      }

      vp = vn;
      ds.params(ds.var_param) = vn(ds.xdim);
      ds.theta = vn(ds.xdim + 1);
    }
    ds.params[ds.inc_param] += ds.delta_inc;
  }
}