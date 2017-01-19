#include <iostream>
#include <vector>
#include <random>
#include "PointProcessFunctions.cpp"

using std::cout;

int main() {
    srand(time(NULL));
    int numberOfSteps;
    cout << "Enter the number of steps in Markov process: ";
    std::cin >> numberOfSteps;
    std::vector<int> markov = generateTwoStepMarkovChain(numberOfSteps);
    cout << "Resulting chain: \n";
    for (int l = 0; l < numberOfSteps - 1; ++l) {
        cout << markov[l] << (l % 30 == 29 ? " \n" : ", ");
    }
    cout << markov[numberOfSteps - 1] << "\n\n";
}



