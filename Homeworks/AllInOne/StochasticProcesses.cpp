//
// Created by tituskc on 11/28/16.
//

#include "StochasticProcesses.h"
#include <iostream>
#include <random>

std::vector<double> StochasticProcesses::generateUnitNormalVector(const int n) const {
    std::vector<double> unitNormals((unsigned long) (n));
    std::mt19937 gen(rand());
    std::normal_distribution<double> norm(0.0, 1.0);
    for (int j = 0; j < n; ++j) {
        unitNormals[j] = norm(gen);
    }
    return unitNormals;
}
