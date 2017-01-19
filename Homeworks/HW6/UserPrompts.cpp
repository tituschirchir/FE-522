//
// Created by tituskc on 11/4/16.
//

#include <iostream>
#include "UserPrompts.h"

using std::cin;
using std::cout;

int UserPrompts::nPrmt() {
    return enterInt("Enter number of steps : ");
}

int UserPrompts::mPrmt() {
    return enterInt("Enter number of simulations: ");
}

int UserPrompts::scalingFactorPrmt() {
    return enterInt("Enter Scaling Factor, K: ");
}

double UserPrompts::timePrmt() {
    return enterDouble("Enter time, T: ");
}

double UserPrompts::timeToCheckPrmt() {
    return enterDouble("Enter time to check. Must be before or at T: ");
}

double UserPrompts::timeFromPrmt() {
    return enterDouble("Enter time from, T_0: ");
}

double UserPrompts::timeToPrmt() {
    return enterDouble("Enter time to, T: ");
}

double UserPrompts::initial() {
    return enterDouble("Enter Initial Value: ");
}

double UserPrompts::mu() {
    return enterDouble("Enter Mean, mu: ");
}

double UserPrompts::sigma() {
    return enterDouble("Enter standard deviation, sigma: ");
}

int UserPrompts::enterInt(const std::string &text) const {
    int N;
    cout << text;
    cin >> N;
    return N;
}

double UserPrompts::enterDouble(const std::string &text) const {
    double N;
    cout << text;
    cin >> N;
    return N;
}

