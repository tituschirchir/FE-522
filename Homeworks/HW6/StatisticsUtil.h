//
// Created by tituskc on 11/3/16.
//

class StatisticalUtil {
    StatisticalUtil();

    ~StatisticalUtil();

public:
    static double variance(const std::vector<double> &values, double mean);

    static double mean(const std::vector<double> &values);
};
