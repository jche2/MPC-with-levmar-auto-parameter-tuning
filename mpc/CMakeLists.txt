cmake_minimum_required(VERSION 3.17)
project(mpc)

include_directories(Eigen-3.3)

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD 11)

add_executable(mpc MPC.cpp cmake-build-debug/mpc_solver.h cmake-build-debug/mpc_solver.cpp MPC.h main.cpp)
target_link_libraries(mpc ipopt ceres)