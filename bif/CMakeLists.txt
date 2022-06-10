project(discrete_bif CXX)
cmake_minimum_required(VERSION 3.21)

find_package (Eigen3 REQUIRED)
find_package (nlohmann_json REQUIRED)

add_executable(main main.cpp dynamical_system.cpp newton.cpp sys_func.cpp)

target_link_libraries (main Eigen3::Eigen nlohmann_json::nlohmann_json)

# target_compile_options(main PUBLIC -O3 -mtune=native -march=native -mfpmath=both)
target_compile_options(main PUBLIC -O3 -mtune=native)