#ifndef CONSTRAINTS_UTILS_H
#define CONSTRAINTS_UTILS_H

typedef double (*funcPtr)(std::vector<double> &x);
typedef bool (*compPtr)(double &x, std::vector<double> &y);

const std::vector<std::string> compSpecial = {"==", ">,<", ">=,<", ">,<=", ">=,<="};
const std::vector<std::string> compHelper = {"<=", "<", "<", "<=", "<="};

// Constraint functions
double prodCpp(std::vector<double> &v);
double sumCpp(std::vector<double> &v);
double meanCpp(std::vector<double> &v);
double maxCpp(std::vector<double> &v);
double minCpp(std::vector<double> &v);

// Standard comparison functions
bool lessCpp(double &x, std::vector<double> &y);
bool greaterCpp(double &x, std::vector<double> &y);
bool lessEqualCpp(double &x, std::vector<double> &y);
bool greaterEqualCpp(double &x, std::vector<double> &y);
bool equalCpp(double &x, std::vector<double> &y);

// In between comparisons
bool greaterLessCpp(double &x, std::vector<double> &y);
bool greaterEqlLessCpp(double &x, std::vector<double> &y);
bool greaterLessEqlCpp(double &x, std::vector<double> &y);
bool greaterEqlLessEqlCpp(double &x, std::vector<double> &y);

Rcpp::XPtr<funcPtr> putFunPtrInXPtr(std::string fstr);
Rcpp::XPtr<compPtr> putCompPtrInXPtr(std::string fstr);

#endif
