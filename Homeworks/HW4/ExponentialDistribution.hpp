
#include <random>
#include <iostream>

class ExponentialDistribution : public Distribution {
    double lambda;
public:
    ExponentialDistribution(unsigned long size, double lambda) : Distribution(size) {
        this->lambda = lambda;
    }

    int getFrequencyToDisplayRatio() { return 10; }

    void populateSample(std::vector<double> &values) {
        std::random_device rd;
        std::mt19937 gen(rd());
        dataPoints["Lambda"] = lambda;
        std::exponential_distribution<> d(lambda);
        for (int n = 0; n < values.size(); ++n) {
            values[n] = d(gen);
        }
    }

    char *getName() {
        return (char *) "Exponential Distribution";
    }

};