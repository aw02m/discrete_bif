#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_

// #include <QString>
#include "qcustomplot.h"
#include <QVector>
#include <eigen3/Eigen/Dense>
#include <nlohmann/json.hpp>
#include <vector>

class dynamical_system {
public:
  dynamical_system(const std::string &json_location);

  Eigen::VectorXd function(const Eigen::VectorXd &x);

  unsigned int xdim;

  // graphic & rk params
  std::vector<unsigned int> axis;
  std::vector<double> xrange;
  std::vector<double> yrange;
  std::vector<double> dparams;
  unsigned int max_plot;

  Eigen::VectorXd x0;
  Eigen::VectorXd p;

  Eigen::VectorXd last_state;
  QVector<QCPGraphData> QCPGsol;

  nlohmann::json json;
  void save_json();
};

#endif