#include "StatisticsUtil.hpp"
#include <random>
#include <iostream>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include <algorithm>

class Distribution {

public:
    unsigned long size;
    std::map<double, double> histogram;
    std::map<std::string, double> dataPoints;

    Distribution(unsigned long size) : size(size) {
    }

    virtual int getFrequencyToDisplayRatio() { return 0; };

    virtual void populateSample(std::vector<double> &vector) {};

    virtual char *getName() { return 0; };

    void processDistribution(std::vector<double> &values, bool normalized, std::ofstream &ofstream, int numberOfBins) {
        this->populateSample(values);

        dataPoints["Sample Mean"] = StatisticsUtil::mean(values);
        dataPoints["Sample Variance"] = StatisticsUtil::variance(values, dataPoints["Sample Mean"]);
        generateHistogram(values, normalized, numberOfBins);
    }

    void printHeader(std::ofstream &ofstream) {
        ofstream << "\n-------------------------------------------------------\n";
        ofstream << getName() << "\n";
        ofstream << "-------------------------------------------------------\n";
        ofstream << "Data Size: " << size << "\n";
        for (auto const &val : dataPoints) {
            ofstream << val.first << " = " << val.second << "\n";
        }
    };


    void generateHistogram(std::vector<double> &values, bool normalized, int numberOfBins) {
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
        for (int i = 0; i < numberOfBins; ++i) {
            double minBin = min + i * binWidth;
            if (normalized) {
                binValues[i] = binValues[i] / (double) vectorSize;
            }
            histogram[minBin] = binValues[i];
        }
    }

    void printHistogram(std::ofstream &ofstream) {
        ofstream << "\nHistogram: \n";
        ofstream << "Min for Bin" << std::setw(30) << "Frequency \n";
        for (auto const &val : this->histogram) {
            unsigned long quant = (unsigned long) (val.second * histogram.size() * getFrequencyToDisplayRatio());
            ofstream << std::setw(12) << val.first << "\t |";
            if (quant > 1)
                ofstream << std::string(quant - 1, '-');
            if (quant > 0)
                ofstream << "* ";
            ofstream << val.second << "\n";

        }
    }
};
