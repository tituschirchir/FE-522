cmake_minimum_required(VERSION 2.6)

project(PointProcesses)

add_definitions(-std=c++11)

add_library(functions SHARED PointProcessFunctions.cpp PointProcessFunctions.hpp StatisticsUtil.hpp)

add_executable(hw5task1 task1.cpp)
add_executable(hw5task2 task2.cpp)
add_executable(hw5task3 task3.cpp)
add_executable(hw5task4 task4.cpp)
add_executable(hw5task5 task5.cpp)