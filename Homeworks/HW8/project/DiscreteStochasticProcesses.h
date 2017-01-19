//
// Created by tituskc on 11/28/16.
//

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <map>
#include "StochasticProcesses.h"
using std::vector;

/*
Markov chain.
Random walk.
 */
class DiscreteStochasticProcesses : public StochasticProcesses {

    std::vector<double> generatePoissonDistribution(double lambda, double T);

    std::vector<double> getInterArrivalTimes(const std::vector<double> arrivalTimes);


    std::vector<int> generateTwoStepMarkovChain(int i);

    std::vector<int> generateThreeStepMarkovChain(int i);

    int count(std::vector<int> array, int value);

    double probability(int n, const int m, bool twoStepMarkov);

    std::vector<int> randomWalk(unsigned long steps);

    std::vector<int> drunkWalk(int n, double probability, std::vector<int> vector1);

    void options();
};
