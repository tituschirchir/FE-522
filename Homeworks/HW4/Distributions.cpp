#include <iostream>
#include <fstream>
#include <map>
#include <random>
#include "Distributions.hpp"
#include "Distribution.hpp"
#include "NormalDistribution.hpp"
#include "UniformDistribution.hpp"
#include "ExponentialDistribution.hpp"

const double PI = std::atan(1.0)*4;

int main() {
    using namespace std;
    ofstream histogramsFile;
    ofstream absoluteErrorFile;
    ofstream dataFile;
    // open files in write mode.
    histogramsFile.open("histograms.txt");
    dataFile.open("data.txt");
    absoluteErrorFile.open("absoluteError.txt");
    absoluteErrorFile << "ABSOLUTE ERROR: \n";
    char choice;
    int dataSize;
    int numberOfBins;
    cout<<"How many data points do you need?: ";
    cin>>dataSize;
    cout<<"How many bins for your histogram (recommended less than sqrt(size))?: ";
    cin>>numberOfBins;
    cout<<"Do you want to execute for N(0,1), U(0,1) and E(1) (Y / N)? ";
    std::cin>>choice;
    if(choice=='Y') {
        doTasks1To3(histogramsFile, dataFile, absoluteErrorFile, (unsigned long) dataSize, numberOfBins);
    }
    else {
        /*Task 4:
         * Request user inputs:
         */
        std::map<double, double> histogram;
        std::cout << "Which distribution? for Normal - N, Exponential - E, Uniform - U? ";
        std::cin >> choice;
        std::vector<double> data((unsigned long) dataSize);
        std::cout << "Do you want a normalized histogram? (Y / N): ";
        char norm;
        std::cin >> norm;
        bool normalized = norm == 'Y';
        switch (choice) {
            case 'N': {
                std::cout << "Enter mean and variance (space delimited): ";
                double arg1, arg2;
                std::cin >> arg1 >> arg2;
                histogram = processNormal(arg1, arg2, normalized, data, histogramsFile, numberOfBins);
                absoluteNormalError(histogram, absoluteErrorFile, dataSize);
                printData(data, dataFile);
                break;
            }
            case 'U': {

                std::cout << "Enter [a, b] (space delimited): ";
                double arg1, arg2;
                std::cin >> arg1 >> arg2;
                histogram = processUniform(arg1, arg2, normalized, data, histogramsFile, numberOfBins);
                absoluteUniformError(histogram, absoluteErrorFile, dataSize);
                printData(data, dataFile);
                break;
            }
            case 'E': {

                std::cout << "Enter Lambda: ";
                double arg1;
                std::cin >> arg1;
                histogram = processExponential(arg1, normalized, data, histogramsFile, numberOfBins);
                absoluteExponentialError(histogram, absoluteErrorFile, dataSize);
                printData(data, dataFile);
                break;
            }
            default: {
                std::cout << choice << " entered is not a valid distribution option. Try E, U or N";
                exit(1);
            }
        }
    }
    // close the opened files.
    histogramsFile.close();
    dataFile.close();
    absoluteErrorFile.close();
    cout<<"\nSee the files: histograms.txt, data.txt and absoluteError.txt for output in root directory\n";
    return 0;
}


void absoluteNormalError(std::map<double, double> values, std::ofstream &ofstream, unsigned long size) {
    double sum = 0.0;
    for (auto const &ent : values) {
        sum+=normalFunction(ent.first);
    }
    double epsilon = 0.0;
    for (auto const &ent : values) {
        double diff = normalFunction(ent.first)/sum - ent.second;
        epsilon += fabs(diff);
    }
    ofstream << "Normal Dist, " << size << " data points: Abs Err: "
              <<epsilon<< " \n";
}

void absoluteExponentialError(std::map<double, double> values, std::ofstream &ofstream, unsigned long size) {
    double sum = 0.0;
    for (auto const &ent : values) {
        sum+=exponentialFunction(ent.first);
    }
    double epsilon = 0.0;
    for (auto const &ent : values) {
        double diff = exponentialFunction(ent.first)/sum - ent.second;
        epsilon += fabs(diff);
    }
    ofstream << "Exponential Dist, " << size << " data points: Abs Err: "
              << epsilon << " \n";
}
void absoluteUniformError(std::map<double, double> values, std::ofstream &ofstream, unsigned long size) {
    double sum = values.size();
    double epsilon = 0.0;
    for (auto const &ent : values) {
        double diff = 1.0/sum - ent.second;
        epsilon += fabs(diff);
    }
    ofstream << "Uniform Dist, " << size << " data points: Abs Err: "
              << epsilon << " \n";
}

void doTasks1To3(std::ofstream &histogramsFile, std::ofstream &dataFile, std::ofstream &absoluteError,
                 unsigned long dataSize, int numberOfBins) {/*Task 1:
     * Generate three vectors of size N (choose 1000) of normally distributed numbers N(0,1), Uniformly distributed
     * random numbers between 0 and 1 U(0,1) and exponentially distributed random numbers E(y) y =1
     * compute the mean and variance of each of the vectors
    */


    /*Task 2
     * Generate Histograms for the vectors generated above. We follow the methodology of binning the generated values in range buckets.
     * For n data points, we choose n^0.5 bins and count the number of items in each bin to get the frequency. e,g for a printNormal distribution which
     * generates 9 numbers 6.5,3,5,4,3.5,6,7,9,8.9,5.5 : # of Bins = 9^0.5 = 3, Bin width= (max-min)/ (# of bins-1) = (9-3)/(3-1) = 2,
     * Bin 1 = {3,4} Freq =  2; normalize = 2/9
     * Bin 2 = {5, 5.5, 6, 6.5,7} Freq = 5; normalize = 5/9
     * Bin 3 = {8.9, 9} Freq = 2; normalize = 2/9
     * Histogram:
     * Min in Bin
     * 3        |-* 0.2222
     * 5        |---* 0.5556
     * 7        |-* 0.2222
     */

    std::vector<double> normalData(dataSize);
    std::vector<double> uniformData(dataSize);
    std::vector<double> exponentialData(dataSize);
    //Task 2: Prints histograms to histograms.txt
    std::map<double, double> normalHistogram = processNormal(0, 1, true, normalData, histogramsFile, numberOfBins);
    std::map<double, double> uniformHistogram = processUniform(0, 1, true, uniformData, histogramsFile, numberOfBins);
    std::map<double, double> exponentialHistogram = processExponential(1, true, exponentialData, histogramsFile,
                                                                       numberOfBins);

    //Task 1: execution: prints data to task1.txt
    printAllData(dataFile, normalData, uniformData, exponentialData);

    /*Task 3:
     * For 10000 data points, we check the first two moments:
     * Normal Distribution, N (0, 1) -> Mean: {-0.00429079 approx 0 as expected}, Variance: {0.977858 approx 1 as expected}
     * Uniform Distribution, U (0, 1) -> Mean: {0.505181 approx 0.5 as expected}, Variance: {0.0839032 approx 0.833333 as expected}
     * Exponential Distribution, E (1) -> Mean: {1.01263 approx 1 as expected}, Variance: {1.05005 approx 1 as expected}
     * The loop below simulates a scenario where we raise size od data by order of ten from 10 to a million
     */
    absoluteNormalError(normalHistogram, absoluteError, dataSize);
    absoluteExponentialError(exponentialHistogram, absoluteError, dataSize);
    absoluteUniformError(uniformHistogram, absoluteError, dataSize);
    /*
     * We observe that in all three distributions, the mean and variance approaches the theoretical value as
     * the size of the sample data increases.
     */
}

void printAllData(std::ofstream &dataFile, const std::vector<double> &normalData,
                  const std::vector<double> &uniformData,
                  const std::vector<double> &exponentialData) {
    dataFile << "Normal Distribution: ";
    printData(normalData, dataFile);

    dataFile << "\nUniform Distribution: ";
    printData(uniformData, dataFile);

    dataFile << "\nExponential Distribution: ";
    printData(exponentialData, dataFile);
}

void printData(std::vector<double> sampleData, std::ofstream &dataFile) {
    double mean = StatisticsUtil::mean(sampleData);
    double variance = StatisticsUtil::variance(sampleData, mean);
    const unsigned long size = sampleData.size();
    dataFile << "Data size = " << size << ", Mean = " << mean << ", Variance = " << variance << "\n";
    for (int i = 0; i < size - 1; i++) {
        if(i % 15==0)
        {
            dataFile<<"\n";
        }
        dataFile << sampleData[i] << ", ";
    }
    dataFile << sampleData[size - 1] << "\n";
}

std::map<double, double>
processExponential(double lambda, bool normalized, std::vector<double> &values, std::ofstream &outfile,
                   int numberOfBins) {
    ExponentialDistribution distribution(values.size(), lambda);
    distribution.processDistribution(values, normalized, outfile, numberOfBins);
    distribution.printHeader(outfile);
    distribution.printHistogram(outfile);
    return distribution.histogram;
}

std::map<double, double> processUniform(double begin, double end, bool normalized, std::vector<double> &values, std::ofstream &outfile,
                                        int numberOfBins) {
    UniformDistribution distribution(values.size(), begin, end);
    distribution.processDistribution(values, normalized, outfile, numberOfBins);
    distribution.printHeader(outfile);
    distribution.printHistogram(outfile);
    return distribution.histogram;
}

std::map<double, double> processNormal(double mean, double variance, bool normalized, std::vector<double> &values, std::ofstream &outfile,
                                       int numberOfBins) {
    NormalDistribution distribution(values.size(), mean, variance);
    distribution.processDistribution(values, normalized, outfile, numberOfBins);
    distribution.printHeader(outfile);
    distribution.printHistogram(outfile);
    return distribution.histogram;
}

double normalFunction(double x) {
    return exp(-1 * x * x / 2) / sqrt(2 * PI);
}

double exponentialFunction(double x) {
    return exp(-1 * x);
}


