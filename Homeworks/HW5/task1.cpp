#include <iostream>
#include <random>
#include <iomanip>
#include <algorithm>
#include <map>
#include "StatisticsUtil.hpp"
#include "PointProcessFunctions.cpp"

using std::cin;
using std::cout;
using std::endl;
using std::vector;


int main() {
    srand(time(NULL));
    cout << "Please enter Lambda and Stopping Time, T: ";
    double lambda, T;
    cin >> lambda >> T;
    vector<double> arrivalTimes = generatePoissonDistribution(lambda, T);
    vector<double> interArrival = getInterArrivalTimes(arrivalTimes);

    double mu = StatisticsUtil::mean(interArrival);
    cout << "InterArrival Times {Sample Mean: " << mu << ", Actual Mean: " << 1 / lambda << "}" << endl;
    double var = StatisticsUtil::variance(interArrival, mu);
    cout << "InterArrival Times {Sample Variance: " << var << ", Actual Variance: " << 1 / (lambda * lambda) << "}"
         << endl;

    printHistogram(generateHistogram(interArrival, true, 20));
    return 0;
}


