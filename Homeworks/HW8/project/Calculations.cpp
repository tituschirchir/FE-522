//
// Created by tituskc on 11/28/16.
//

#include <algorithm>
#include <iomanip>
#include "Calculations.h"


double Calculations::max(double one, double two) {
    return one > two ? one : two;
}

double Calculations::d2(double S, double K, double T, double r, double vol) {
    return (log(S / K) + (r - vol * vol) * T) / (vol * sqrt(T));
}

void Calculations::separator() { printf("---------------------------\n"); }

double Calculations::mean(const std::vector<double> &values) {
    double total = 0.0;
    for (int i = 0; i < values.size(); i++) {
        total += values[i];
    }
    return total / values.size();
}
double Calculations::variance(const std::vector<double> &values, double mean) {
    double totalVariation = 0.0;
    for (int i = 0; i < values.size(); ++i) {
        double variationFromMean = values[i] - mean;
        totalVariation += variationFromMean * variationFromMean;
    }
    return totalVariation / (values.size() - 1);
}

std::map<double, double> Calculations::buildHistogram(std::vector<double> &values, bool normalized, int numberOfBins) {
    std::map<double, double> histogram;
    std::sort(values.begin(), values.end());
    const unsigned long vectorSize = values.size();
    double max = values[vectorSize - 1], min = values[0];

    const double diff = max - min;
    double binWidth = diff / numberOfBins;

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
    return histogram;
}
void Calculations::printHistogram(std::ofstream &ofstream, std::map<double, double> histogram){
    ofstream << "\nHistogram: \n";
    ofstream << "Min for Bin" << std::setw(30) << "Frequency \n";
    for (auto const &val : histogram) {
        unsigned long quant = (unsigned long) (val.second * histogram.size());
        ofstream << std::setw(12) << val.first << "\t |";
        if (quant > 1)
            ofstream << std::string(quant - 1, '-');
        if (quant > 0)
            ofstream << "* ";
        ofstream << val.second << "\n";
    }
}

double Calculations::randomWithDefaultZeroToOne(double min, double max) {
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

void Calculations::writeToFile(const std::vector<double> &data, const std::string fileName, const double dt) {
    std::ofstream myfile;
    myfile.open(fileName + ".xls");
    for (int i = 0; i < data.size(); ++i) {
        myfile << i * dt << "," << data[i] << "\n";
    }
    myfile.close();

    std::printf("Output: See %s.xls at the base of directory.\n", fileName.c_str());
}