//
// Created by tituskc on 11/3/16.
//
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "BrownianMotion.h"
#include "StatisticsUtil.h"

using std::vector;

BrownianMotion::BrownianMotion() {}

BrownianMotion::~BrownianMotion() {}

vector<double> BrownianMotion::scaledRandomWalk(const double from, const double to, const int k) {
    if (to < from) {
        std::cout << from << " - " << to << " is an invalid range. Please select from less than to" << std::endl;
        exit(1);
    }
    double tk = k * (to - from);
    double scalingRoot = 1 / sqrt(k);
    vector<double> randomWalk((unsigned long) tk);
    int delta_i = 0;
    for (int i = 1; i < tk; ++i) {
        delta_i = (double) rand() / RAND_MAX <= 0.5 ? -1 : 1;
        randomWalk[i] = scalingRoot * (randomWalk[i - 1] + delta_i);
    }
    return randomWalk;
}

vector<double> BrownianMotion::standardBrownianMotion(double T, int N, double initial) {
    vector<double> myVec((unsigned long) N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> norm(0, 1);
    myVec[0] = initial;
    double dt = T / N;
    for (int i = 1; i < N; i++) {
        myVec[i] = myVec[i - 1] + sqrt(dt) * norm(gen);
    }
    return myVec;
}

vector<double>
BrownianMotion::brownianMotionWithDrift(const double mu, const double sigma, const double T, const int n,
                                        const double initial) {
    double dt = T / n;
    vector<double> myVec((unsigned long) n);
    vector<double> unitNormals = generateUnitNormalVector(n);
    myVec[0] = initial;
    for (int i = 1; i < n; i++) {
        myVec[i] = myVec[i - 1] + (sigma * sqrt(dt) * unitNormals[i] + mu * dt);
    }
    return myVec;
}

vector<double>
BrownianMotion::geometricBrownianMotion(const double S0, const double mu, const double sigma, const double T,
                                        const int N) {
    double dt = T / N;
    vector<double> geoBrown(N);
    vector<double> unitNormals = generateUnitNormalVector(N);
    geoBrown[0] = S0;
    for (int i = 0; i < N-1; i++) {
        geoBrown[i+1] = geoBrown[i] * exp(sigma * sqrt(dt) * unitNormals[i+1] + mu * dt);
    }
    return geoBrown;
}

void BrownianMotion::checkMeanAndVarianceWithDrift(const double sigma, const double mu, const double T,
                                                   const double timeToCheck, const double N, int M) {
    if (timeToCheck > T) {
        std::cout << "Please select a time less than or equal to T";
        exit(1);
    }
    vector<vector<double>> bm = multidimensionalBrownian(mu, sigma, T, N, M);
    vector<double> timeSlice(M);
    for (int i = 0; i < M; ++i) {
        timeSlice[i] = bm[i][timeToCheck * N / T - 1];
    }
    double mean = StatisticalUtil::mean(timeSlice);
    double variance = StatisticalUtil::variance(timeSlice, mean);
    std::cout << "\n---------------------------------------------------------------------------" << std::endl;
    std::cout << "B(t+dt) - B(t) is normally distributed with mean mu*t and Variance sigma^2*T\n";
    std::cout << "Theoretical E[B(t)] = " << (T * mu) << "; Empirical E[B(t)] = " << mean << "\n";
    std::cout << "Theoretical Var[B(t)] = " << (T * sigma * sigma) << "; Empirical Var[B(t)] = " << variance << "\n";
    std::cout << "---------------------------------------------------------------------------" << std::endl;
}

void BrownianMotion::checkMeanAndVarianceWithoutDrift(const double T, const double timeToCheck, const double N, int M) {
    if (timeToCheck > T) {
        std::cout << "Please select a time less than or equal to T";
        exit(1);
    }
    vector<vector<double>> bm = multidimensionalStandardBrownian(T, N, M);
    vector<double> timeSlice((unsigned long) M);
    for (int i = 0; i < M; ++i) {
        timeSlice[i] = bm[i][timeToCheck * N / T - 1];
    }
    double mean = StatisticalUtil::mean(timeSlice);
    double variance = StatisticalUtil::variance(timeSlice, mean);
    std::cout << "\n---------------------------------------------------------------------------" << std::endl;
    std::cout << "Theoretical E[B(t)] = " << 0 << "; Empirical E[B(t)] = " << mean << "\n";
    std::cout << "Theoretical Var[B(t)] = " << timeToCheck << "; Empirical Var[B(t)] = " << variance << "\n";
    std::cout << "---------------------------------------------------------------------------" << std::endl;
}

vector<vector<double>> BrownianMotion::multidimensionalStandardBrownian(double T, int N, int M) {
    vector<vector<double>> grid(M, vector<double>(N));
    for (int i = 0; i < M; ++i) {
        grid[i] = standardBrownianMotion(T, N, 0);
    }
    return grid;
}

vector<vector<double>> BrownianMotion::multidimensionalBrownian(double mu, double sigma, double T, int N, int M) {
    vector<vector<double>> grid(M, vector<double>(N));
    for (int i = 0; i < M; ++i) {
        grid[i] = brownianMotionWithDrift(mu, sigma, T, N, 0);
    }
    return grid;
}

std::vector<std::vector<double>> BrownianMotion::twoDBrownianMotion(double mu, double sigma, double T, int N) {
    return (std::vector<std::vector<double>>) {brownianMotionWithDrift(mu, sigma, T, N, 0),
                                               brownianMotionWithDrift(mu, sigma, T, N, 0)};
}

void BrownianMotion::writeArrayToFile(const vector<vector<double>> &array, const std::string &fileName) {
    std::ofstream myfile;
    myfile.open(fileName);
    for (int i = 0; i < array[0].size(); ++i) {
        myfile << i;
        for (int j = 0; j < array.size(); ++j) {
            myfile << "," << array[j][i];
        }
        myfile << "\n";
    }
    myfile.close();
}

void BrownianMotion::writeArrayToFile(const vector<double> &array, const std::string &fileName) {
    std::ofstream myfile;
    myfile.open(fileName);
    for (int i = 0; i < array.size(); ++i) {
        myfile << i << "," << array[i] << "\n";
    }
    myfile.close();
}

vector<double> BrownianMotion::generateUnitNormalVector(const int n) const {
    vector<double> unitNormals((unsigned long) (n));
    std::mt19937 gen(rand());
    std::normal_distribution<double> norm(0.0, 1.0);
    for (int j = 0; j < n; ++j) {
        unitNormals[j] = norm(gen);
    }
    return unitNormals;
}

