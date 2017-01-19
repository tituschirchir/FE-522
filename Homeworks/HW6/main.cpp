#include <iostream>
#include <vector>
#include "BrownianMotion.h"
#include "UserPrompts.h"
#include "homework6.h"

BrownianMotion bm;
UserPrompts up;
using std::vector;

/*General Requirements
 * To construct and handle arrays, you can use either pointers (to
 * pointers for multidimensional arrays) or the class vector. Build your own class, and name it
 * BrownianMotion. It should be regarded as a simulator for processes related with BM. Put
 * function declarations into .hpp files in your project directories and include your .hpp files in
 * your .cpp files. Put your function definitions into separate .cpp files so that the one with your
 * main function only has function calls. Also, please include in your main file all the .cpp files
 * used for definitions.
 */

int main() {
    std::cout << "Which task do you want to execute? For task 1 (Qn 1 in Homework) enter 1, task 2: 2 etc... : ";
    int task = 0;
    std::cin >> task;
    switch (task) {
        case 1:
            task1();
            break;
        case 2:
            task2();
            break;
        case 3:
            task3();
            break;
        case 4:
            task4();
            break;
        case 5:
            task5();
            break;
        case 6:
            task6();
            break;
        case 7:
            task7();
            break;
        default:
            std::cout << "Incorrect input: Enter a number between 1 and 7.\n";
            exit(1);
    }
    return 0;
}

void task1() {
    std::cout
            << "Welcome to the Scaled Random Walk. For simulations, you will be required to enter 'From Time, T0', 'To Time, T'"
                    " and the scaling factor, k. We recommend a large k and require T>T0: Enjoy!! \n";
    const char *fileName = "ascaledrandomwalk.xls";
    bm.writeArrayToFile(bm.scaledRandomWalk(up.timeFromPrmt(), up.timeToPrmt(), up.scalingFactorPrmt()), fileName);
    std::cout << "Done! Output file can be found at ./" << fileName;
}

void task2() {
    std::cout
            << "Simulation of Standard Brownian Motion simulation. You can enter inputs at your convenience. These include"
                    " Time T, Number of steps N and initial value B0. Enjoy.\n\n";
    const char *fileName = "astandardbrownianmotion.xls";
    bm.writeArrayToFile(bm.standardBrownianMotion(up.timePrmt(), up.nPrmt(), up.initial()), fileName);
    std::cout << "Done! Output file can be found at ./" << fileName;

}

void task3() {
    std::cout
            << "Simulation of Brownian Motion with Drift. Enter Time T, Number of steps N and initial value B0. Enjoy.\n\n";
    const char *fileName = "abrownianmotionwithdrift.xls";
    bm.writeArrayToFile(bm.brownianMotionWithDrift(up.mu(), up.sigma(), up.timePrmt(), up.nPrmt(), up.initial()),
                        fileName);
    std::cout << "Done! Output file can be found at ./" << fileName;
}

void task4() {
    std::cout
            << "Welcome to Geometric Brownian Motion simulation. You can enter inputs at your convenience. "
                    "These include Time T, Number of steps N and initial value B0. Enjoy.\n\n";
    const char *fileName = "ageometricbrownianmotion.xls";
    bm.writeArrayToFile(bm.geometricBrownianMotion(up.initial(), up.mu(), up.sigma(), up.timePrmt(), up.nPrmt()),
                        fileName);
    std::cout << "Done! Output file can be found at ./" << fileName;
}

void task5() {
    std::cout << "For this task, we will test the properties of Brownian motion without drift. That is, E[B(t)] = 0"
            " and Var[B(t)] = t. You will be required to enter T, time at which to check, number of steps N and Number of simulations, M\n\n";
    bm.checkMeanAndVarianceWithoutDrift(up.timePrmt(), up.timeToCheckPrmt(), up.nPrmt(), up.mPrmt());
}

void task6() {
    std::cout << "For this task, we will test the properties of Brownian motion without drift. That is, E[B(t)] = μt"
            " and Var[B(t)] = σ^2 t. You will be required to enter T, time at which to check, number of steps N and Number of simulations, M\n\n";
    bm.checkMeanAndVarianceWithDrift(up.sigma(), up.mu(), up.timePrmt(), up.timeToCheckPrmt(), up.nPrmt(), up.mPrmt());
}

void task7() {
    std::cout
            << "For this task, we will generate 2 D Brownian motion (essentially two brownian motion paths) see hw6Report.pdf "
                    "for the plots. You will be prompted to enter Mean mu, Volatility, sigma, Time T and Number of steps, N. Enjoy.\n\n";
    const char *fileName = "atwoDbrownianmotion.xls";
    bm.writeArrayToFile(bm.twoDBrownianMotion(up.mu(), up.sigma(), up.timePrmt(), up.nPrmt()), fileName);
    std::cout << "Done! Output file can be found at ./" << fileName;
}

