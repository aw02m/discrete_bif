#include "dynamical_system.hpp"
#include <fstream>
#include <iostream>

dynamical_system::dynamical_system(const std::string &json_location) {
  std::ifstream ifs(json_location);
  if (ifs.fail()) {
    std::cerr << "File does NOT exist." << std::endl;
    std::exit(1);
  }
  ifs >> json;

  xdim = json["x0"].size();

  axis = json["axis"].get<std::vector<unsigned int>>();
  xrange = json["xrange"].get<std::vector<double>>();
  yrange = json["yrange"].get<std::vector<double>>();
  dparams = json["dparams"].get<std::vector<double>>();
  max_plot = json["max_plot"];

  out_path = json["out_path"];
  json_out_path = json["json_out_path"];

  /* These json array should be casted to the STL container type*/
  std::vector<double> fixed_arr = json["x0"];
  Eigen::Map<Eigen::VectorXd> x0(fixed_arr.data(), fixed_arr.size());
  this->x0 = x0;
  this->last_state = x0;

  std::vector<double> params_arr = json["params"];
  Eigen::Map<Eigen::VectorXd> params(params_arr.data(), params_arr.size());
  this->p = params;

  trajectory_set =
      std::vector<Eigen::VectorXd>(max_plot, Eigen::VectorXd::Zero(xdim));
}

void dynamical_system::save_json() {
  json["x0"] = last_state;
  json["params"] = p;
  std::ofstream json_out;
  json_out.open(json_out_path, std::ios::out);
  json_out << json.dump(4);
  json_out.close();
  Eigen::IOFormat Out(Eigen::FullPrecision, 0, " ", "\n", " ", " ");
  std::ofstream f;
  f.open(out_path, std::ios::out);
  f << std::fixed;
  for (unsigned int i = 0; i < max_plot; i++) {
    f << trajectory_set[i].transpose().format(Out) << std::endl;
  }
  f.close();
}