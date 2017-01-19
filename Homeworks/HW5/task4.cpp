#include <iostream>
#include <vector>
#include <random>
#include "PointProcessFunctions.cpp"

using std::cout;
using std::endl;
int main() {
    srand(time(NULL));
    int numberOfSteps;
    cout << "Enter the number of steps in Markov process: ";
    std::cin >> numberOfSteps;
    std::vector<int> markov = generateThreeStepMarkovChain(numberOfSteps);
    cout << "Markov chain: \n";
    for (int l = 0; l < numberOfSteps - 1; ++l) {
        cout << markov[l] << (l % 30 == 29 ? " \n" : ", ");
    }
    cout << markov[numberOfSteps - 1] << "\n\n";
    cout << "The probability of getting state 1 for 1,000,000 steps is " << probability(10000000, 1, false)<<endl;
}



