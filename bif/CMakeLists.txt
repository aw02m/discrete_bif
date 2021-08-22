project(discrete_bif CXX)
include_directories("C:\\Program Files (x86)\\Eigen3\\include")
add_executable(main main.cpp dynamical_system.cpp ds_derivatives.cpp ds_func.cpp newton.cpp)