cmake_minimum_required(VERSION 3.17)

project(discrete_bif LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### INCLUSION OF Numerical Libs
find_package (Eigen3 REQUIRED)
find_package (nlohmann_json REQUIRED)

### INCLUSION OF Project Files
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")

include_directories(${INCLUDE_DIR})
include_directories(${SOURCE_DIR})

file(GLOB_RECURSE SOURCES
    "${INCLUDE_DIR}/*.h"
    "${SOURCE_DIR}/*.cpp"
)

### SETUP EXECUTABLE
add_executable(main ${SOURCES})

target_include_directories(main PRIVATE ${INCLUDE_DIR})
target_include_directories(main PRIVATE ${SOURCE_DIR})

target_link_libraries (main PRIVATE
  Eigen3::Eigen
  nlohmann_json::nlohmann_json
)

# target_compile_options(main PUBLIC -O3 -mtune=native -march=native -mfpmath=both)
target_compile_options(main PUBLIC -O3 -mtune=native)