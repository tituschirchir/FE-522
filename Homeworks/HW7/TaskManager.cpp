//
// Created by tituskc on 11/12/16.
//
#include "TaskManager.h"

#include <random>
#include <fstream>
#include <boost/math/distributions/normal.hpp>

const std::string VASICEK = "Vasicek";
const std::string CIR = "CIR";
using std::vector;

TaskManager::TaskManager() {}

TaskManager::~TaskManager() {}

std::vector<double>
TaskManager::euroMonteCarloSimulation(double S, double mu, double sigma, double T, int N, int M, double K) {
    std::vector<double> mcSimulation(M);
    for (int i = 0; i < M; ++i) {
        mcSimulation[i] = max(geometricBrownianMotion(S, mu, sigma, T, N)[N - 1] - K, 0);
    }
    return mcSimulation;
}

std::vector<double>
TaskManager::asianMonteCarloSimulation(double S, double mu, double sigma, double T, int N, int M, double K) {
    vector<double> mcSimulation(M);
    vector<double> values;
    for (int i = 0; i < M; ++i) {
        values = geometricBrownianMotion(S, mu, sigma, T, N);
        mcSimulation[i] = max(mean(values) - K, 0);
    }
    return mcSimulation;
}

std::vector<double>
TaskManager::simulateIRModel(const std::string model, double a, double b, double r0, double vol, double T, int N) {
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

double TaskManager::europeanCallOptionPriceBlackScholes(double S, double K, double T, double r, double vol) {
    double c = d2(S, K, T, r, vol);
    double nd1 = cumulativeNormalDistribution(c + vol * sqrt(T));
    double nd2 = cumulativeNormalDistribution(c);
    return S * nd1 - exp(-1 * r * T) * K * nd2;
}

double
TaskManager::europeanCallOptionPriceMonteCarlo(double S, double sigma, double T, int N, int M, double K, double r) {
    double mu = r - sigma * sigma / 2;
    vector<double> mcSimulation = euroMonteCarloSimulation(S, mu, sigma, T, N, M, K);
    double price = exp(-1 * r * T) * mean(mcSimulation);
    return price;
}

double TaskManager::asianCallOptionPriceMonteCarlo(double S, double sigma, double T, int N, int M, double K, double r) {
    double mu = r - sigma * sigma / 2;
    vector<double> mcSimulation = asianMonteCarloSimulation(S, mu, sigma, T, N, M, K);
    return exp(-1 * r * T) * mean(mcSimulation);
}

double
TaskManager::zeroCouponBondPrice(const std::string model, double a, double b, double r0, double vol, double T, int N,
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
void TaskManager::execute() {
    int mission = 0;
    int displayOptions = 1;
    do {
        missionControl(mission, displayOptions);
        separator();
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
        separator();
    } while (mission != 0);

}

void TaskManager::missionControl(int &mission, int &displayOptions) const {
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

void TaskManager::options() const {
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

void TaskManager::pricingEuropeanCallOptionUsingBlackScholes() {

    std::printf("GOAL: Price European Call using Black Scholes\n");
    double S, K, r, vol, T;
    std::cout << "Please enter S, K, T, r, vol. (space delimited - order matters): ";
    std::cin >> S >> K >> T >> r >> vol;
    std::printf("Inputs: S=%g, K=%g, T=%g, r=%g, vol=%g)\n", S, K, T, r, vol);
    double euroCallPrice = europeanCallOptionPriceBlackScholes(S, K, T, r, vol);
    std::printf("Price: %g \n", euroCallPrice);

}

void TaskManager::pricingEuropeanCallOptionUsingMonteCarloSimulation() {

    std::printf("GOAL: Price European Call using Monte Carlo Simulation\n");
    double S, K, r, vol, T;
    int N, M;
    std::cout << "Please enter S, K, T, r, vol, N and M (space delimited - order matters): ";
    std::cin >> S >> K >> T >> r >> vol >> N >> M;
    std::printf("Inputs: S=%g, K=%g, T=%g, r=%g , vol=%g, N=%d, M=%d \n", S, K, T, r, vol, N, M);
    double euroCallPrice = europeanCallOptionPriceMonteCarlo(S, vol, T, N, M, K, r);
    std::printf("Price: %g \n", euroCallPrice);
}

void TaskManager::checkEuroCallPricingByMCvsBS() {

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

void TaskManager::pricingAsianCallByMonteCarloSimulation() {

    std::cout << "GOAL: Price Asian Option using Monte Carlo Simulation\n";
    double S, K, r, vol, T;
    int N, M;
    std::cout << "Please enter S, K, T, r, vol, N and M (space delimited - order matters): ";
    std::cin >> S >> K >> T >> r >> vol >> N >> M;
    std::printf("Inputs : S=%g, K=%g, T=%g, r=%g , vol=%g, N=%d, M=%d \n", S, K, T, r, vol, N, M);
    double price = asianCallOptionPriceMonteCarlo(S, vol, T, N, M, K, r);
    std::printf("Price %g\n", price);
}

void TaskManager::generateOrnsteinUhlenbeckInterestRateProcess() {

    std::printf("GOAL: Generate Ornstein-Uhlenbeck interest rate process\n");
    double a, b, r0, vol, T;
    int N;
    std::cout << "Please enter a, b, r0, vol, T and N (space delimited - order matters): ";
    std::cin >> a >> b >> r0 >> vol >> T >> N;
    std::printf("Inputs :: a=%g, b=%g, r0=%g, vol=%g, T=%g, N=%d \n", a, b, r0, vol, T, N);
    write(simulateIRModel("Vasicek", a, b, r0, vol, T, N), "vasicek", T / N);
}

void TaskManager::generateCoxIngersolRossCoxInterestRateProcess() {

    std::printf("GOAL: Generate Cox-Ingersol-Ross interest rate process\n");
    double a, b, r0, vol, T;
    int N;
    std::cout << "Please enter a, b, r0, vol, T and N (space delimited - order matters): ";
    std::cin >> a >> b >> r0 >> vol >> T >> N;
    std::printf("Inputs :: a=%g, b=%g, r0=%g, vol=%g, T=%g, N=%d \n", a, b, r0, vol, T, N);
    write(simulateIRModel("CIR", a, b, r0, vol, T, N), "cir", T / N);
}

void TaskManager::pricingZeroCouponBondUsingMonteCarloSimulation() {

    std::printf("GOAL: Price Zero Coupon Bond bu Monte Carlo simulation using CIR interest rate model\n");
    double a, b, r0, vol, T;
    int N;
    int M;
    std::cout << "Please enter a, b, r0, vol, T, N and M(Monte carlo simulations) (space delimited - order matters): ";
    std::cin >> a >> b >> r0 >> vol >> T >> N >> M;
    std::printf("Inputs: a=%g, b=%g, r0=%g, vol=%g, T=%g, N=%d, M=%d \n", a, b, r0, vol, T, N, M);
    std::printf("Price: %g \n", zeroCouponBondPrice("CIR", a, b, r0, vol, T, N, M));
}

void TaskManager::separator() const { printf("---------------------------\n"); }

double TaskManager::mean(const std::vector<double> values) {
    double total = 0.0;
    for (int i = 0; i < values.size(); i++) {
        total += values[i];
    }
    return total / values.size();
}

double TaskManager::cumulativeNormalDistribution(double x) {
    boost::math::normal norm;
    return 1 - cdf(complement(norm, x));
}

double TaskManager::max(double one, double two) {
    return one > two ? one : two;
}

double TaskManager::d2(double S, double K, double T, double r, double vol) {
    return (log(S / K) + (r - vol * vol) * T) / (vol * sqrt(T));
}

vector<double> TaskManager::scaledUnitNormalVector(const int n, const double deltaT) const {
    vector<double> unitNormals((unsigned long) (n));
    std::mt19937 gen(rand());
    std::normal_distribution<double> norm(0.0, 1.0);
    for (int j = 0; j < n; ++j) {
        unitNormals[j] = norm(gen) * sqrt(deltaT);
    }
    return unitNormals;
}

vector<double> TaskManager::geometricBrownianMotion(const double S0, const double mu, const double sigma,
                                                    const double T, const int N) {
    vector<double> geoBrown = scaledUnitNormalVector(N, T / N);
    geoBrown[0] = S0;
    for (int i = 1; i < N; i++) {
        geoBrown[i] = geoBrown[i - 1] * exp(sigma * geoBrown[i] + mu * (T / N));
    }
    return geoBrown;
}

void TaskManager::write(const std::vector<double> &data, const std::string fileName, const double dt) const {
    std::ofstream myfile;
    myfile.open(fileName + ".xls");
    for (int i = 0; i < data.size(); ++i) {
        myfile << i * dt << "," << data[i] << "\n";
    }
    myfile.close();

    std::printf("Output: See %s.xls at the base of directory.\n", fileName.c_str());
}
