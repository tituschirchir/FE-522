
#include <random>
#include <iostream>

class UniformDistribution : public Distribution {
    double begin;
    double end;
public:
    UniformDistribution(unsigned long size, double begin, double end) : Distribution(size) {
        this->begin = begin;
        this->end = end;
    }

    int getFrequencyToDisplayRatio() { return 50; }

    void populateSample(std::vector<double> &values) {
        std::random_device rd;
        std::mt19937 gen(rd());
        dataPoints["Begin"] = begin;
        dataPoints["End"] = end;

        std::uniform_real_distribution<> dis(begin, end);
        for (int n = 0; n < values.size(); ++n) {
            values[n] = dis(gen);
        }
    }

    char *getName() {
        return (char *) "Uniform Distribution";
    }

};