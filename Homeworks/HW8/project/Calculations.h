//
// Created by tituskc on 11/28/16.
//
#include <vector>
#include <fstream>
#include <map>

class Calculations {
public:
    static double variance(const std::vector<double> &values, double mean);

    static double mean(const std::vector<double> &values);

    std::map<double, double> buildHistogram(std::vector<double> &values, bool normalized, int numberOfBins);

    void printHistogram(std::ofstream &ofstream, std::map<double, double> histogram);

    static void writeToFile(const std::vector<double> &data, const std::string fileName, const double dt);

    static double d2(double S, double K, double T, double r, double vol);

    static double max(double one, double two);

    static void separator();

    static double randomWithDefaultZeroToOne(double min = 0.0, double max = 1.0);
};