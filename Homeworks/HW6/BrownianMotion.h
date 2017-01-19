//
// Created by tituskc on 11/3/16.
//
#include <vector>

class BrownianMotion {

public:
    BrownianMotion();

    ~BrownianMotion();

    /*1.
     * In the previous assignment, you had to simulate a well–known Markov Chain, namely the
     * Random Walk. Generate now a scaled simple symmetric Random Walk {B_k (t) : t ≥ 0} on
     * an interval [0, T ] chosen by the user. Note that the scaling factor k ∈ N ∗ should be a large
     * integer. Hint. Please see the Appendix 1.5 in the attached file, BM.pdf.
     */
    std::vector<double> scaledRandomWalk(const double from, const double to, const int k);

    /*2.
     * A direct application of the Central Limit Theorem (CLT) yields B_k (t) → B(t) ∼ N (0, t) as
     * k → ∞, in distribution. Simulate a standard BM at times 0 = t_0 < t_1 < t_2 < · · · < t_N = T ,
     * using samples from a standard normal random variable. Simulate also a BM starting at
     * T, i = 0, 1, . . . , N. With some a ∈ R, i.e. B(0) = a. Indications. Pick for simplicity t_i = i*T/N
     * this, we have ∆t = t i+1 −t i = N , ∀i. Your program should output a csv file with two vectors:
     * the time vector and the simulated BM. Use Microsoft Excel to plot your results.
     */
    std::vector<double> standardBrownianMotion(double T, int N, double initial);

    /*3.
     * Simulate a BM with drift μ ∈ R and variance term σ > 0, X(t) = σB(t) + μt, on an interval
     * [0, T ] chosen by the user. Use the points 0 = t 0 < t 1 < t 2 < · · · < t_N = T , and output a
     * csv file containing the time vector and X(t k ), k = 0, 1, . . . , N , as before. Plot your results
     * in Microsoft Excel.
     */
    std::vector<double>
    brownianMotionWithDrift(const double mu, const double sigma, const double T, const int n, const double initial);

    /*4.
     * Simulate now a Geometric BM S(t) = S(0)e^X(t) , t ≥ 0, where X(t) = σB(t) + μt.
     * Use Microsoft Excel to plot your results.
     */
    std::vector<double>
    geometricBrownianMotion(const double S0, const double mu, const double sigma, const double T, const int N);

    /*5.
     * Now, we shall check some theoretical properties of the BM and related processes. Generate
     * M = 1 million paths of a BM, denoted B_1 (t), B_2 (t), . . . , B_M (t), with t ∈ [0, T ]. Check
     * numerically that E[B(T)] = 0, Var(B(T)) = T. Obviously, the same relation must hold for each
     * t_k , i.e. E[B(t_k )] = 0, Var(B(t k )) = t k . Hint. The main tool we use is a Monte Carlo Simulation.
     */

    void checkMeanAndVarianceWithoutDrift(const double T, const double timeToCheck, const double N, int M);

    /*6.
    * Check numerically that the BM with drift X(t) = σB(t) + μt has expectation and variance
    * given by E[X(t)] = μt, Var(X(t)) = σ^2 t.
    */
    void checkMeanAndVarianceWithDrift(const double sigma, const double mu, const double T,
                                       const double timeToCheck, const double N, int M);

    /* 7.
     * Simulate and plot a 2 Dimensional Brownian Motion. Hint. Generate two BMs, one for each dimension,
     * say B_x (t) and B_y (t). Output the vectors into a csv file and simply plot(B_x , B_y )
     */
    std::vector<std::vector<double>> twoDBrownianMotion(double mu, double sigma, double T, int N);

    void writeArrayToFile(const std::vector<std::vector<double>> &array, const std::string &fileName);

    void writeArrayToFile(const std::vector<double> &array, const std::string &fileName);

private:
    std::vector<std::vector<double>> multidimensionalBrownian(double mu, double sigma, double T, int N, int M);

    std::vector<std::vector<double>> multidimensionalStandardBrownian(double T, int N, int M);

    std::vector<double> generateUnitNormalVector(const int n) const;
};
