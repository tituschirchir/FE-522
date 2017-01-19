#include <iostream>
#include <random>
#include <map>
#include <algorithm>
#include <iomanip>
#include "PointProcessFunctions.hpp"

using std::cin;
using std::cout;
using std::endl;
using std::vector;
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

vector<double> getInterArrivalTimes(const vector<double> arrivalTimes) {
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

vector<double> generatePoissonDistribution(double lambda, double T) {
    vector<double> arrivalTimes;
    double t = 0.0;
    while (t < T) {
        double nextTime = (-1 / lambda) * log(randomWithDefaultZeroToOne());
        t = t + nextTime;
        arrivalTimes.push_back(t);
    }
    return arrivalTimes;
}

map<double, double> generateHistogram(vector<double> &values, bool normalized, int numberOfBins) {
    std::sort(values.begin(), values.end());
    const unsigned long vectorSize = values.size();
    double max = values[vectorSize - 1];
    double min = values[0];

    const double diff = max - min;
    double binWidth = diff / (numberOfBins);

    double *binValues = new double[numberOfBins];
    for (int i = 0; i < vectorSize; i++) {
        ++binValues[(int) ((values[i] - min) / binWidth)];
    }
    map<double, double> histogram;
    for (int i = 0; i < numberOfBins; ++i) {
        double minBin = min + i * binWidth;
        if (normalized) {
            binValues[i] = binValues[i] / (double) vectorSize;
        }
        histogram[minBin] = binValues[i];
    }
    return histogram;
}

void printHistogram(const map<double, double> histogram) {
    cout << "\nHistogram: \n";
    cout << "Min for Bin" << std::setw(30) << "Frequency \n";
    double multiplier = histogram.size() * 8 / histogram.at(0);
    for (auto const &val : histogram) {
        unsigned long quant = (unsigned long) (val.second * multiplier);
        cout << std::setw(12) << val.first << "\t |";
        if (quant > 1)
            cout << std::string(quant - 1, '-');
        cout << "* " << val.second << "\n";

    }
}


std::vector<int> generateTwoStepMarkovChain(int numberOfSteps) {
    int X = S2D[0];
    std::vector<int> markov;
    for (int k = 0; k < numberOfSteps; ++k) {
        markov.push_back(X);
        double random = randomWithDefaultZeroToOne();
        if (X == S2D[0] && random > P[0][0]) {
            X = S2D[1];
        } else if (X == S2D[1] && random <= P[1][0]) {
            X = S2D[0];
        }
    }
    return markov;
}

std::vector<int> generateThreeStepMarkovChain(int numberOfSteps) {
    int X = S3D[0];
    std::vector<int> markov;
    markov.push_back(X);
    for (int k = 0; k < numberOfSteps; ++k) {
        double random = randomWithDefaultZeroToOne();
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

double probability(int n, const int m, bool twoStepMarkov) {
    double sum = 0.0;
    for (int i = 0; i < m; i++) {
        std::vector<int> markov = twoStepMarkov ? generateTwoStepMarkovChain(n) : generateThreeStepMarkovChain(n);
        int num = count(markov, 1);
        sum += (double) num / n;
    }
    return sum / m;
}

int count(std::vector<int> array, int value) {
    int c = 0;
    for (int i = 0; i < array.size(); ++i) {
        if (array[i] == value)
            c++;
    }
    return c;
}

double randomWithDefaultZeroToOne(double min, double max) {
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}


vector<int> randomWalk(unsigned long steps) {
    vector<int> v(steps);
    return drunkWalk(0, 0.666666666666666, v);
}

vector<int> drunkWalk(int n, double probability, vector<int> vector1) {
    if (n >= (int) vector1.size()) {
        return vector1;
    } else {
        double type = randomWithDefaultZeroToOne();
        int prev = n == 0 ? 0 : vector1[n - 1];
        vector1[n] = type < probability ? prev + 1 : prev - 1;
        return drunkWalk(++n, probability, vector1);
    }
}
