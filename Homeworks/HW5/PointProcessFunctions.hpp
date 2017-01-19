
#include <map>
#include <vector>

void printHistogram(const std::map<double, double> hist);

std::map<double, double> generateHistogram(std::vector<double> &values, bool normalized, int numberOfBins);

std::vector<double> generatePoissonDistribution(double lambda, double T);

std::vector<double> getInterArrivalTimes(const std::vector<double> arrivalTimes);


std::vector<int> generateTwoStepMarkovChain(int i);

std::vector<int> generateThreeStepMarkovChain(int i);

int count(std::vector<int> array, int value);

double probability(int n, const int m, bool twoStepMarkov);

double randomWithDefaultZeroToOne(double min = 0.0, double max = 1.0);

std::vector<int> randomWalk(unsigned long steps);

std::vector<int> drunkWalk(int n, double probability, std::vector<int> vector1);


