#include "dynamical_system.hpp"
#include "newton.hpp"
#include <nlohmann/json.hpp>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Put a path to the input json." << std::endl;
    std::exit(1);
  }

  std::ifstream ifs(argv[1]);
  if (ifs.fail()) {
    std::cerr << "File does NOT exist." << std::endl;
    std::exit(1);
  }

  nlohmann::json json;
  ifs >> json;
  dynamical_system ds(json);
  if (ds.mode < 0 || ds.mode > 3) {
    std::cerr << "Specified mode does not exist." << std::endl;
    std::cerr << "Set \"mode\" in JSON as FIX:0, G:1, PD:2, NS:3" << std::endl;
    std::exit(1);
  }

  std::cout << "**************************************************"
            << std::endl;
  std::cout << "<initial condition>" << std::endl;
  std::cout << "mode : ";
  switch (ds.mode) {
  case 0:
    std::cout << "FIXED POINT" << std::endl;
    break;
  case 1:
    std::cout << "TANGENT BIF" << std::endl;
    break;
  case 2:
    std::cout << "PERIOD DOUBLING BIF" << std::endl;
    break;
  case 3:
    std::cout << "NEIMARK-SACKER BIF" << std::endl;
    break;
  }
  std::cout << "increment parameter : " << ds.inc_param << std::endl;
  std::cout << "system dimention : " << ds.xdim << std::endl;
  std::cout << "period : " << ds.period << std::endl;
  std::cout << "params  : ";
  std::cout << ds.p.transpose() << std::endl;
  std::cout << "x0  : ";
  std::cout << ds.x0.transpose() << std::endl;
  std::cout << "**************************************************"
            << std::endl;

  newton(ds);

  // output json for latest state
  json["x0"] = ds.x0;
  json["params"] = ds.p;
  std::ofstream json_out;
  json_out.open(ds.json_out_path, std::ios::out);
  json_out << json.dump(4);
  json_out.close();

  return 0;
}