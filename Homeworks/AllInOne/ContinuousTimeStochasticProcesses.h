//
// Created by tituskc on 11/28/16.
//

#include <vector>
#include <iostream>
#include "StochasticProcesses.h"

class ContinuousTimeStochasticProcesses : public StochasticProcesses {
    // poisson, ohrenstein-uhnbeck, brownian motion
    std::vector<double> standardBrownianMotion(double T, int N, double initial);

    std::vector<double>
    brownianMotionWithDrift(const double mu, const double sigma, const double T, const int n, const double initial);

    std::vector<double> geometricBrownianMotion(const double S0, const double mu, const double sigma, const double T, const int N);

    void checkMeanAndVarianceWithoutDrift(const double T, const double timeToCheck, const double N, int M);
    void checkMeanAndVarianceWithDrift(const double sigma, const double mu, const double T,
                                       const double timeToCheck, const double N, int M);

    std::vector<std::vector<double>> twoDBrownianMotion(double mu, double sigma, double T, int N);

    void writeArrayToFile(const std::vector<std::vector<double>> &array, const std::string &fileName);

    void writeArrayToFile(const std::vector<double> &array, const std::string &fileName);

    void missionControl(int &mission, int &displayOptions);

    double europeanCallOptionPriceBlackScholes(double forwardSpotPrice, double strikePrice, double yearsToExpiry,
                                               double rate, double vol);

    double europeanCallOptionPriceMonteCarlo(double S, double sigma, double T, int N, int M, double K, double r);

    double asianCallOptionPriceMonteCarlo(double S, double sigma, double T, int N, int M, double K, double r);

    double
    zeroCouponBondPrice(const std::string model, double a, double b, double r0, double vol, double T, int N, int M);

    double cumulativeNormalDistribution(double x);

    std::vector<double> euroMonteCarloSimulation(double S, double mu, double sigma, double T, int N, int M, double K);

    std::vector<double> asianMonteCarloSimulation(double S, double mu, double sigma, double T, int N, int M, double K);

    std::vector<double> simulateIRModel(const std::string model, double a, double b, double r0, double vol, double T, int N);

    std::vector<double> scaledUnitNormalVector(const int n, const double deltaT) const;

    void pricingEuropeanCallOptionUsingBlackScholes();

    void pricingEuropeanCallOptionUsingMonteCarloSimulation();

    void checkEuroCallPricingByMCvsBS();

    void pricingAsianCallByMonteCarloSimulation();

    void generateOrnsteinUhlenbeckInterestRateProcess();

    void generateCoxIngersolRossCoxInterestRateProcess();

    void pricingZeroCouponBondUsingMonteCarloSimulation();

    std::vector<std::vector<double>> multidimensionalBrownian(double mu, double sigma, double T, int N, int M);

    std::vector<std::vector<double>> multidimensionalStandardBrownian(double T, int N, int M);

public:
    void execute();

    void options();
};
