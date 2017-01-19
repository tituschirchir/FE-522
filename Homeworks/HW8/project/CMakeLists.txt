cmake_minimum_required(VERSION 2.6)

project(Distributions)

add_definitions(-std=c++11)

add_executable(stochasticprocesses main.cpp
        StochasticProcesses.cpp StochasticProcesses.h
        DiscreteStochasticProcesses.cpp DiscreteStochasticProcesses.h
        ContinuousTimeStochasticProcesses.cpp ContinuousTimeStochasticProcesses.h
        Calculations.cpp Calculations.h)
