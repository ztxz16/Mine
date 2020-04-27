cmake_minimum_required(VERSION 3.5)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -pthread --std=c++11")

add_executable(Solver Solver.cpp)
add_executable(MineSolver MineSolver.cpp)