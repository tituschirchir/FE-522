cmake_minimum_required(VERSION 2.6)
project(PricingProcesses)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp TaskManager.cpp TaskManager.h)
add_executable(h7 ${SOURCE_FILES})
