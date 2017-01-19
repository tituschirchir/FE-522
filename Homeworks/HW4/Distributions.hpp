std::map<double, double>
processNormal(double mean, double variance, bool normalized, std::vector<double> &values, std::ofstream &outfile,
              int numberOfBins);

std::map<double, double>
processUniform(double begin, double end, bool normalized, std::vector<double> &values, std::ofstream &outfile,
               int numberOfBins);

std::map<double, double>
processExponential(double lambda, bool normalized, std::vector<double> &values, std::ofstream &outfile,
                   int numberOfBins);

double normalFunction(double x);

double exponentialFunction(double x);

void printData(std::vector<double> sampleData, std::ofstream &dataFile);

void printAllData(std::ofstream &dataFile, const std::vector<double> &normalData,
                  const std::vector<double> &uniformData,
                  const std::vector<double> &exponentialData);

void doTasks1To3(std::ofstream &histogramsFile, std::ofstream &dataFile, std::ofstream &ofstream, unsigned long dataSize,
                 int numberOfBins);

void absoluteUniformError(std::map<double, double> values, std::ofstream &ofstream, unsigned long i);

void absoluteExponentialError(std::map<double, double> values, std::ofstream &ofstream, unsigned long i);

void absoluteNormalError(std::map<double, double> values, std::ofstream &ofstream, unsigned long i);
