cmake_minimum_required(VERSION 2.6)

project(BrownianMotion)

add_definitions(-std=c++11)

add_executable(hw6 main.cpp homework6.h BrownianMotion.cpp BrownianMotion.h UserPrompts.cpp UserPrompts.h StatisticsUtil.h StatisticsUtil.cpp)
