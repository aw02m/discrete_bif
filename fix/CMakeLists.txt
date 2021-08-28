project(discrete_fix CXX)
include_directories("C:\\Program Files (x86)\\Eigen3\\include")
add_executable(main main.cpp dynamical_system.cpp eigensolver.cpp newton.cpp ds_func.cpp ds_derivatives.cpp)