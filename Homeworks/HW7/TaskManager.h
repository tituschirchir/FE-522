//
// Created by tituskc on 11/12/16.
//

#ifndef PRICINGPROCESSES_PRICING_H
#define PRICINGPROCESSES_PRICING_H


#include <iostream>
#include <vector>

class TaskManager {
public:
    TaskManager();

    ~TaskManager();

    void execute();

private:

    void options() const;

    void missionControl(int &mission, int &displayOptions) const;

    double europeanCallOptionPriceBlackScholes(double forwardSpotPrice, double strikePrice, double yearsToExpiry,
                                               double rate, double vol);

    double europeanCallOptionPriceMonteCarlo(double S, double sigma, double T, int N, int M, double K, double r);

    double asianCallOptionPriceMonteCarlo(double S, double sigma, double T, int N, int M, double K, double r);

    double
    zeroCouponBondPrice(const std::string model, double a, double b, double r0, double vol, double T, int N, int M);

    double mean(const std::vector<double> values);

    double max(double one, double two);

    double cumulativeNormalDistribution(double x);

    std::vector<double> euroMonteCarloSimulation(double S, double mu, double sigma, double T, int N, int M, double K);

    std::vector<double> asianMonteCarloSimulation(double S, double mu, double sigma, double T, int N, int M, double K);

    std::vector<double>
    simulateIRModel(const std::string model, double a, double b, double r0, double vol, double T, int N);

    std::vector<double> scaledUnitNormalVector(const int n, const double deltaT) const;

    std::vector<double>
    geometricBrownianMotion(const double S0, const double mu, const double sigma, const double T, const int N);

    double d2(double forwardSpotPrice, double strikePrice, double yearsToExpiry, double rate, double vol);

    void pricingEuropeanCallOptionUsingBlackScholes();

    void pricingEuropeanCallOptionUsingMonteCarloSimulation();

    void checkEuroCallPricingByMCvsBS();

    void pricingAsianCallByMonteCarloSimulation();

    void generateOrnsteinUhlenbeckInterestRateProcess();

    void generateCoxIngersolRossCoxInterestRateProcess();

    void pricingZeroCouponBondUsingMonteCarloSimulation();

    void write(const std::vector<double> &vac, const std::string fileName, const double dt) const;

    void separator() const;
};


#endif //PRICINGPROCESSES_PRICING_H
