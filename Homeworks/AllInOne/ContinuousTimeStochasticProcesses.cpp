//
// Created by tituskc on 11/28/16.
//

#include "ContinuousTimeStochasticProcesses.h"
#include "Calculations.h"
#include <cmath>
#include <random>
#include <boost/math/distributions/normal.hpp>
using namespace std;
const std::string VASICEK = "Vasicek";
const std::string CIR = "CIR";

vector<double> ContinuousTimeStochasticProcesses::standardBrownianMotion(double T, int N, double initial) {
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
ContinuousTimeStochasticProcesses::brownianMotionWithDrift(const double mu, const double sigma, const double T, const int n,
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
ContinuousTimeStochasticProcesses::geometricBrownianMotion(const double S0, const double mu, const double sigma, const double T,
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

void ContinuousTimeStochasticProcesses::checkMeanAndVarianceWithDrift(const double sigma, const double mu, const double T,
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
    double mean = Calculations::mean(timeSlice);
    double variance = Calculations::variance(timeSlice, mean);
    std::cout << "\n---------------------------------------------------------------------------" << std::endl;
    std::cout << "B(t+dt) - B(t) is normally distributed with mean mu*t and Variance sigma^2*T\n";
    std::cout << "Theoretical E[B(t)] = " << (T * mu) << "; Empirical E[B(t)] = " << mean << "\n";
    std::cout << "Theoretical Var[B(t)] = " << (T * sigma * sigma) << "; Empirical Var[B(t)] = " << variance << "\n";
    std::cout << "---------------------------------------------------------------------------" << std::endl;
}

void ContinuousTimeStochasticProcesses::checkMeanAndVarianceWithoutDrift(const double T, const double timeToCheck, const double N, int M) {
    if (timeToCheck > T) {
        std::cout << "Please select a time less than or equal to T";
        exit(1);
    }
    vector<vector<double>> bm = multidimensionalStandardBrownian(T, N, M);
    vector<double> timeSlice((unsigned long) M);
    for (int i = 0; i < M; ++i) {
        timeSlice[i] = bm[i][timeToCheck * N / T - 1];
    }
    double mean = Calculations::mean(timeSlice);
    double variance = Calculations::variance(timeSlice, mean);
    std::cout << "\n---------------------------------------------------------------------------" << std::endl;
    std::cout << "Theoretical E[B(t)] = " << 0 << "; Empirical E[B(t)] = " << mean << "\n";
    std::cout << "Theoretical Var[B(t)] = " << timeToCheck << "; Empirical Var[B(t)] = " << variance << "\n";
    std::cout << "---------------------------------------------------------------------------" << std::endl;
}

vector<vector<double>> ContinuousTimeStochasticProcesses::multidimensionalStandardBrownian(double T, int N, int M) {
    vector<vector<double>> grid(M, vector<double>(N));
    for (int i = 0; i < M; ++i) {
        grid[i] = standardBrownianMotion(T, N, 0);
    }
    return grid;
}

vector<vector<double>> ContinuousTimeStochasticProcesses::multidimensionalBrownian(double mu, double sigma, double T, int N, int M) {
    vector<vector<double>> grid(M, vector<double>(N));
    for (int i = 0; i < M; ++i) {
        grid[i] = brownianMotionWithDrift(mu, sigma, T, N, 0);
    }
    return grid;
}

std::vector<std::vector<double>> ContinuousTimeStochasticProcesses::twoDBrownianMotion(double mu, double sigma, double T, int N) {
    return (std::vector<std::vector<double>>) {brownianMotionWithDrift(mu, sigma, T, N, 0),
                                               brownianMotionWithDrift(mu, sigma, T, N, 0)};
}

void ContinuousTimeStochasticProcesses::writeArrayToFile(const vector<vector<double>> &array, const std::string &fileName) {
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

void ContinuousTimeStochasticProcesses::writeArrayToFile(const vector<double> &array, const std::string &fileName) {
    std::ofstream myfile;
    myfile.open(fileName);
    for (int i = 0; i < array.size(); ++i) {
        myfile << i << "," << array[i] << "\n";
    }
    myfile.close();
}


std::vector<double>
ContinuousTimeStochasticProcesses::euroMonteCarloSimulation(double S, double mu, double sigma, double T, int N, int M, double K) {
    std::vector<double> mcSimulation(M);
    for (int i = 0; i < M; ++i) {
        mcSimulation[i] = Calculations::max(geometricBrownianMotion(S, mu, sigma, T, N)[N - 1] - K, 0);
    }
    return mcSimulation;
}

std::vector<double>
ContinuousTimeStochasticProcesses::asianMonteCarloSimulation(double S, double mu, double sigma, double T, int N, int M, double K) {
    vector<double> mcSimulation(M);
    vector<double> values;
    for (int i = 0; i < M; ++i) {
        values = geometricBrownianMotion(S, mu, sigma, T, N);
        mcSimulation[i] = max(Calculations::mean(values) - K, 0.0);
    }
    return mcSimulation;
}

std::vector<double>
ContinuousTimeStochasticProcesses::simulateIRModel(const std::string model, double a, double b, double r0, double vol, double T, int N) {
    if (VASICEK != model && CIR != model) {
        std::cout << "Please select " << VASICEK << " or " << CIR << " model. " << std::endl;
        exit(1);
    }
    vector<double> vec = scaledUnitNormalVector(N, T / N);
    vector<double> cir(N);
    cir[0] = r0;
    for (int i = 1; i <= N; ++i) {
        cir[i] = cir[i - 1] + a * (b - cir[i - 1]) * T / N + vol * (CIR == model ? sqrt(cir[i - 1]) : 1) * vec[i];
    }
    return cir;
}

double ContinuousTimeStochasticProcesses::europeanCallOptionPriceBlackScholes(double S, double K, double T, double r, double vol) {
    double c = Calculations::d2(S, K, T, r, vol);
    double nd1 = cumulativeNormalDistribution(c + vol * sqrt(T));
    double nd2 = cumulativeNormalDistribution(c);
    return S * nd1 - exp(-1 * r * T) * K * nd2;
}

double
ContinuousTimeStochasticProcesses::europeanCallOptionPriceMonteCarlo(double S, double sigma, double T, int N, int M, double K, double r) {
    double mu = r - sigma * sigma / 2;
    vector<double> mcSimulation = euroMonteCarloSimulation(S, mu, sigma, T, N, M, K);
    double price = exp(-1 * r * T) * Calculations::mean(mcSimulation);
    return price;
}

double ContinuousTimeStochasticProcesses::asianCallOptionPriceMonteCarlo(double S, double sigma, double T, int N, int M, double K, double r) {
    double mu = r - sigma * sigma / 2;
    vector<double> mcSimulation = asianMonteCarloSimulation(S, mu, sigma, T, N, M, K);
    return exp(-1 * r * T) * Calculations::mean(mcSimulation);
}

double
ContinuousTimeStochasticProcesses::zeroCouponBondPrice(const std::string model, double a, double b, double r0, double vol, double T, int N,
                                 int M) {
    vector<double> interestRates;
    double sum = 0.0;
    double reimannSum;
    for (int j = 0; j < M; ++j) {
        reimannSum = 0.0;
        interestRates = simulateIRModel(model, a, b, r0, vol, T, N);
        for (int i = 0; i < N; ++i) {
            reimannSum += interestRates[i] * T / N;
        }
        sum += exp(-1 * reimannSum);
    }
    return sum / M;
}

//User prompts
void ContinuousTimeStochasticProcesses::execute() {
    int mission = 0;
    int displayOptions = 1;
    do {
        missionControl(mission, displayOptions);
        Calculations::separator();
        switch (mission) {
            case 0:
                std::cout << "Good Bye! The total processing time for your exploits was ";
                continue;
            case 1:
                pricingEuropeanCallOptionUsingBlackScholes();
                break;
            case 2:
                pricingEuropeanCallOptionUsingMonteCarloSimulation();
                break;
            case 3:
                checkEuroCallPricingByMCvsBS();
                break;
            case 4:
                pricingAsianCallByMonteCarloSimulation();
                break;
            case 5:
                generateOrnsteinUhlenbeckInterestRateProcess();
                break;
            case 6:
                generateCoxIngersolRossCoxInterestRateProcess();
                break;
            case 7:
                pricingZeroCouponBondUsingMonteCarloSimulation();
                break;
            default:
                break;
        }
        Calculations::separator();
    } while (mission != 0);

}

void ContinuousTimeStochasticProcesses::missionControl(int &mission, int &displayOptions) {
    if (mission == 0) {
        options();
        displayOptions = 0;
    } else {
        std::cout << "To see the options again, enter 1. Otherwise enter 0: ";
        std::cin >> displayOptions;
        if (displayOptions == 1) {
            options();
            displayOptions = 0;
        }
    }
    std::cout << "Choose Mission: ";
    std::cin >> mission;
}

void ContinuousTimeStochasticProcesses::options() {
    std::cout << "Select Mission: \n"
            "   1: Price European Call using Black-Scholes\n"
            "   2: Price European Call using Monte Carlo Simulation\n"
            "   3: Compare Monte Carlo and Black-Scholes Option Prices\n"
            "   4: Price Asian Call using Monte Carlo Simulation\n"
            "   5: Generate data dump of interest rate process using Ornstein-Uhlenbeck model\n"
            "   6: Generate data dump of interest rate process using Cox-Ingersol-Ross(CIR) model\n"
            "   7: Price Zero coupon bond using Monte Carlo and CIR for interest rate process\n"
            "   0: Exit\n";
}

void ContinuousTimeStochasticProcesses::pricingEuropeanCallOptionUsingBlackScholes() {

    std::printf("GOAL: Price European Call using Black Scholes\n");
    double S, K, r, vol, T;
    std::cout << "Please enter S, K, T, r, vol. (space delimited - order matters): ";
    std::cin >> S >> K >> T >> r >> vol;
    std::printf("Inputs: S=%g, K=%g, T=%g, r=%g, vol=%g)\n", S, K, T, r, vol);
    double euroCallPrice = europeanCallOptionPriceBlackScholes(S, K, T, r, vol);
    std::printf("Price: %g \n", euroCallPrice);

}

void ContinuousTimeStochasticProcesses::pricingEuropeanCallOptionUsingMonteCarloSimulation() {

    std::printf("GOAL: Price European Call using Monte Carlo Simulation\n");
    double S, K, r, vol, T;
    int N, M;
    std::cout << "Please enter S, K, T, r, vol, N and M (space delimited - order matters): ";
    std::cin >> S >> K >> T >> r >> vol >> N >> M;
    std::printf("Inputs: S=%g, K=%g, T=%g, r=%g , vol=%g, N=%d, M=%d \n", S, K, T, r, vol, N, M);
    double euroCallPrice = europeanCallOptionPriceMonteCarlo(S, vol, T, N, M, K, r);
    std::printf("Price: %g \n", euroCallPrice);
}

void ContinuousTimeStochasticProcesses::checkEuroCallPricingByMCvsBS() {

    std::printf("GOAL: Compare Black-Scholes and Monte Carlo Simulation in Pricing European Call Option\n");
    double S, K, r, vol, T;
    int N, M;
    std::cout << "Please enter S, K, T, r, vol, N and M (space delimited - order matters): ";
    std::cin >> S >> K >> T >> r >> vol >> N >> M;
    std::printf("Inputs : S=%g, K=%g, T=%g, r=%g , vol=%g, N=%d, M=%d \n", S, K, T, r, vol, N, M);
    double bsPrice = europeanCallOptionPriceBlackScholes(S, K, T, r, vol);
    std::printf("Black-Scholes Price = %g \n", bsPrice);

    double mcPrice = europeanCallOptionPriceMonteCarlo(S, vol, T, N, M, K, r);
    std::printf("Monte Carlo Price = %g \n", mcPrice);

    double diff = fabs(mcPrice - bsPrice);
    std::cout << "Difference: " << diff << " = " << diff * 100 / bsPrice<<"%"<< std::endl;

}

void ContinuousTimeStochasticProcesses::pricingAsianCallByMonteCarloSimulation() {

    std::cout << "GOAL: Price Asian Option using Monte Carlo Simulation\n";
    double S, K, r, vol, T;
    int N, M;
    std::cout << "Please enter S, K, T, r, vol, N and M (space delimited - order matters): ";
    std::cin >> S >> K >> T >> r >> vol >> N >> M;
    std::printf("Inputs : S=%g, K=%g, T=%g, r=%g , vol=%g, N=%d, M=%d \n", S, K, T, r, vol, N, M);
    double price = asianCallOptionPriceMonteCarlo(S, vol, T, N, M, K, r);
    std::printf("Price %g\n", price);
}

void ContinuousTimeStochasticProcesses::generateOrnsteinUhlenbeckInterestRateProcess() {

    std::printf("GOAL: Generate Ornstein-Uhlenbeck interest rate process\n");
    double a, b, r0, vol, T;
    int N;
    std::cout << "Please enter a, b, r0, vol, T and N (space delimited - order matters): ";
    std::cin >> a >> b >> r0 >> vol >> T >> N;
    std::printf("Inputs :: a=%g, b=%g, r0=%g, vol=%g, T=%g, N=%d \n", a, b, r0, vol, T, N);
    Calculations::writeToFile(simulateIRModel("Vasicek", a, b, r0, vol, T, N), "vasicek", T / N);
}

void ContinuousTimeStochasticProcesses::generateCoxIngersolRossCoxInterestRateProcess() {

    std::printf("GOAL: Generate Cox-Ingersol-Ross interest rate process\n");
    double a, b, r0, vol, T;
    int N;
    std::cout << "Please enter a, b, r0, vol, T and N (space delimited - order matters): ";
    std::cin >> a >> b >> r0 >> vol >> T >> N;
    std::printf("Inputs :: a=%g, b=%g, r0=%g, vol=%g, T=%g, N=%d \n", a, b, r0, vol, T, N);
    Calculations::writeToFile(simulateIRModel("CIR", a, b, r0, vol, T, N), "cir", T / N);
}

void ContinuousTimeStochasticProcesses::pricingZeroCouponBondUsingMonteCarloSimulation() {

    std::printf("GOAL: Price Zero Coupon Bond bu Monte Carlo simulation using CIR interest rate model\n");
    double a, b, r0, vol, T;
    int N;
    int M;
    std::cout << "Please enter a, b, r0, vol, T, N and M(Monte carlo simulations) (space delimited - order matters): ";
    std::cin >> a >> b >> r0 >> vol >> T >> N >> M;
    std::printf("Inputs: a=%g, b=%g, r0=%g, vol=%g, T=%g, N=%d, M=%d \n", a, b, r0, vol, T, N, M);
    std::printf("Price: %g \n", zeroCouponBondPrice("CIR", a, b, r0, vol, T, N, M));
}


double ContinuousTimeStochasticProcesses::cumulativeNormalDistribution(double x) {
    boost::math::normal norm;
    return 1 - cdf(complement(norm, x));
}

vector<double> ContinuousTimeStochasticProcesses::scaledUnitNormalVector(const int n, const double deltaT) const {
    vector<double> unitNormals((unsigned long) (n));
    std::mt19937 gen(rand());
    std::normal_distribution<double> norm(0.0, 1.0);
    for (int j = 0; j < n; ++j) {
        unitNormals[j] = norm(gen) * sqrt(deltaT);
    }
    return unitNormals;
}
