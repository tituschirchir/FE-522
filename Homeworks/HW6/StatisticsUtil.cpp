//
// Created by tituskc on 11/3/16.
//

#include <vector>
#include "StatisticsUtil.h"

StatisticalUtil::StatisticalUtil() {};

StatisticalUtil::~StatisticalUtil() {}

double StatisticalUtil::mean(const std::vector<double> &values) {
    double total = 0.0;
    for (int i = 0; i < values.size(); i++) {
        total += values[i];
    }
    return total / values.size();
}

double StatisticalUtil::variance(const std::vector<double> &values, double mean) {

    double totalVariation = 0.0;
    for (int i = 0; i < values.size(); ++i) {
        double variationFromMean = values[i] - mean;
        totalVariation += variationFromMean * variationFromMean;
    }
    return totalVariation / (values.size() - 1);

}

