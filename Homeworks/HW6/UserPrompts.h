//
// Created by tituskc on 11/4/16.
//

class UserPrompts {
public:
    int nPrmt();

    int mPrmt();

    double timePrmt();

    int scalingFactorPrmt();

    double initial();

    double timeToPrmt();

    double timeFromPrmt();

    double timeToCheckPrmt();

    double mu();

    double sigma();

private:

    int enterInt(const std::string &text) const;

    double enterDouble(const std::string &text) const;
};

