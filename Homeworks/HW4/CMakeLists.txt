cmake_minimum_required(VERSION 2.6)

project(Distributions)

add_definitions(-std=c++11)

add_library(functions SHARED Distributions.hpp Distributions.cpp Distribution.hpp runDistribution.sh StatisticsUtil.hpp ExponentialDistribution.hpp NormalDistribution.hpp UniformDistribution.hpp)

add_executable(hw4 Distributions.cpp)