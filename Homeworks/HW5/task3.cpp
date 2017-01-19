#include "PointProcessFunctions.cpp"

int main() {
    srand(time(NULL));
    cout << "The frequency of appearance of state 1 for  a two state Markov chain for different step sizes: \n";
    cout << "M" << std::setw(10) << "N " << std::setw(20) << "Probability" << endl;
    for (int i = 0; i < 7; i++) {
        int n = (int) powf(10, i);
        cout << n + 1 << std::setw(10) << n << std::setw(20) << probability(n, true, n + 1) << endl;
    }
}
