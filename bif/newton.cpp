#include "newton.hpp"
#include "dynamical_system.hpp"
#include <fstream>
#include <iomanip>
#include <vector>

void newton(dynamical_system &ds) {
  unsigned int target_dim;
  if (ds.mode != 0) {
    target_dim = ds.xdim + 1;
  } else {
    target_dim = ds.xdim;
  }
  Eigen::VectorXd vp(target_dim);
  Eigen::VectorXd vn(target_dim);
  Eigen::VectorXd F(target_dim);
  Eigen::MatrixXd J(target_dim, target_dim);
  double norm;
  Eigen::IOFormat Comma(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]");
  std::cout << std::fixed << std::setprecision(16);
  bool exit_flag = false;

  vp(Eigen::seqN(0, ds.xdim)) = ds.x0;
  if (ds.mode != 0) {
    vp(ds.xdim) = ds.p(ds.var_param);
  }

  Eigen::IOFormat Out(Eigen::FullPrecision, 0, " ", "\n", " ", " ");
  std::ofstream f;
  f.open(ds.out_path, std::ios::out);
  f << std::fixed;

  for (int p = 0; p < ds.inc_iter; p++) {
    if (exit_flag) {
      ds.inc_iter = p - 1;
      break;
    }
    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < ds.max_iter; i++) {
      auto FJ = ds.newton_FJ(vp);
      F = std::get<0>(FJ);
      J = std::get<1>(FJ);
      vn = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>(J).solve(-F) + vp;
      // debug(F);
      // exit_flag = true;

      norm = F.norm();
      if (norm < ds.eps) {
        auto end = std::chrono::system_clock::now();
        auto dur = end - start;
        auto msec =
            std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
        // std::cout << "**************************************************"
        //           << std::endl;
        std::cout << std::setprecision(16);
        std::cout << p << " : converged (iter = " << i + 1 << ", ";
        std::cout << "time = " << msec << "[msec])" << std::endl;
        std::cout << "params : " << ds.p.transpose().format(Comma) << std::endl;
        std::cout << "x0     : "
                  << vn(Eigen::seqN(0, ds.xdim)).transpose().format(Comma)
                  << std::endl;
        std::cout << std::setprecision(4);
        std::cout << "(Re(μ), Im(μ)), abs(μ), arg(μ) :" << std::endl;
        for (int k = 0; k < ds.xdim; k++) {
          std::cout << ds.eigvals(k) << ", ";
          std::cout << std::abs(ds.eigvals(k)) << ", ";
          std::cout << std::arg(ds.eigvals(k)) << std::endl;
        }
        std::cout << "**************************************************"
                  << std::endl;
        f << ds.p.transpose().format(Out) << std::endl;
        vp = vn;
        break;
      } else if (norm >= ds.explode) {
        std::cerr << "explode (iter = " << i + 1 << ")" << std::endl;
        exit_flag = true;
        break;
      }
      if (i >= ds.max_iter - 1) {
        std::cerr << "iter over (iter = " << i + 1 << ")" << std::endl;
        exit_flag = true;
        break;
      }
      vp = vn;
    }
    ds.p[ds.inc_param] += ds.delta_inc;
  }

  // set last state
  ds.x0 = vn(Eigen::seqN(0, ds.xdim));
  ds.p(ds.var_param) = vn(ds.xdim);
  f.close();
}