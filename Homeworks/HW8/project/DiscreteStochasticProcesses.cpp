//
// Created by tituskc on 11/28/16.
//

#include <algorithm>
#include <iomanip>
#include "DiscreteStochasticProcesses.h"
#include "Calculations.h"

using std::cout;
using std::map;
double P[2][2] =
        {
                {0.3, 0.7},
                {0.5, 0.5}
        };
int S2D[2] = {0, 1};

double P3D[3][3] =
        {
                {0.3333333,     0.3333333,     0.3333334},
                {0.125,         0.5,           0.375},
                {0.42857142857, 0.57142857143, 0}
        };

int S3D[3] = {0, 1, 2};

vector<double> DiscreteStochasticProcesses::getInterArrivalTimes(const vector<double> arrivalTimes) {
    unsigned long n = arrivalTimes.size();
    cout << "We have " << n << " stopping times \n";
    vector<double> interArrival(n);
    for (int i = 0; i < n - 2; i++) {
        if (i != n - 2) {
            interArrival[i + 1] = arrivalTimes[i + 1] - arrivalTimes[i];
        }
    }
    return interArrival;
}

vector<double> DiscreteStochasticProcesses::generatePoissonDistribution(double lambda, double T) {
    vector<double> arrivalTimes;
    double t = 0.0;
    while (t < T) {
        double nextTime = (-1 / lambda) * log(Calculations::randomWithDefaultZeroToOne());
        t = t + nextTime;
        arrivalTimes.push_back(t);
    }
    return arrivalTimes;
}


std::vector<int> DiscreteStochasticProcesses::generateTwoStepMarkovChain(int numberOfSteps) {
    int X = S2D[0];
    std::vector<int> markov;
    for (int k = 0; k < numberOfSteps; ++k) {
        markov.push_back(X);
        double random = Calculations::randomWithDefaultZeroToOne();
        if (X == S2D[0] && random > P[0][0]) {
            X = S2D[1];
        } else if (X == S2D[1] && random <= P[1][0]) {
            X = S2D[0];
        }
    }
    return markov;
}

std::vector<int> DiscreteStochasticProcesses::generateThreeStepMarkovChain(int numberOfSteps) {
    int X = S3D[0];
    std::vector<int> markov;
    markov.push_back(X);
    for (int k = 0; k < numberOfSteps; ++k) {
        double random = Calculations::randomWithDefaultZeroToOne();
        if (X == S3D[0] && random > P3D[0][0]) {
            X = S3D[random <= (P3D[0][0] + P3D[0][1]) ? 1 : 2];
        } else if (X == S3D[1]) {
            if (random <= P3D[1][0]) {
                X = S3D[0];
            } else if (random > (P3D[1][0] + P3D[1][1])) {
                X = S3D[2];
            }
        } else if (random <= (P3D[2][0] + P3D[2][1])) {
            X = S3D[(random <= P3D[2][0] ? 0 : 1)];
        }
        markov.push_back(X);
    }
    return markov;
}

double DiscreteStochasticProcesses::probability(int n, const int m, bool twoStepMarkov) {
    double sum = 0.0;
    for (int i = 0; i < m; i++) {
        std::vector<int> markov = twoStepMarkov ? generateTwoStepMarkovChain(n) : generateThreeStepMarkovChain(n);
        int num = count(markov, 1);
        sum += (double) num / n;
    }
    return sum / m;
}

int DiscreteStochasticProcesses::count(std::vector<int> array, int value) {
    int c = 0;
    for (int i = 0; i < array.size(); ++i) {
        if (array[i] == value)
            c++;
    }
    return c;
}


vector<int> DiscreteStochasticProcesses::randomWalk(unsigned long steps) {
    vector<int> v(steps);
    return drunkWalk(0, 0.666666666666666, v);
}
void DiscreteStochasticProcesses::options() {
    cout<<"You have selected the Discrete processes.";
}
vector<int> DiscreteStochasticProcesses::drunkWalk(int n, double probability, vector<int> vector1) {
    if (n >= (int) vector1.size()) {
        return vector1;
    } else {
        double type = Calculations::randomWithDefaultZeroToOne();
        int prev = n == 0 ? 0 : vector1[n - 1];
        vector1[n] = type < probability ? prev + 1 : prev - 1;
        return drunkWalk(++n, probability, vector1);
    }
}