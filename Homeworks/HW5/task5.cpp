
#include <random>
#include <iostream>
#include "PointProcessFunctions.cpp"

using std::cin;
using std::cout;
using std::endl;
using std::vector;

int main() {
    unsigned long numberOfSteps;
    srand(time(NULL));
    cout << "Enter the number of steps for random walk: ";
    cin >> numberOfSteps;
    vector<int> randomVals = randomWalk(numberOfSteps);
    cout << "\nTHE WALK:\n{ ";
    for (int i = 0; i < numberOfSteps - 1; ++i) {
        cout << randomVals[i];
        cout << (i % 30 == 29 ? " \n" : ", ");
    }
    cout << randomVals[numberOfSteps - 1] << " }\n\n";
    return 0;
}





