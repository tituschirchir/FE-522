
#include <random>
#include <iostream>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>

class NormalDistribution : public Distribution {
    double mean;
    double variance;
public:
    NormalDistribution(unsigned long size, double mean, double variance) : Distribution(size) {
        this->mean = mean;
        this->variance = variance;
    }

    int getFrequencyToDisplayRatio() { return 30; }

    void populateSample(std::vector<double> &values) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> d(mean, sqrt(variance));
        dataPoints["Mean"] = mean;
        dataPoints["Variance"] = variance;
        for (int n = 0; n < values.size(); ++n) {
            values[n] = d(gen);
        }
    }

    char *getName() {
        return (char *) "Normal Distribution";
    }

};