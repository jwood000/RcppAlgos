#ifndef CONSTRAINTS_UTILS_H
#define CONSTRAINTS_UTILS_H

typedef double (*funcPtr)(std::vector<double>& x);
typedef bool (*compPtr)(double& x, double& y);

// Constraint functions
double prodCpp(std::vector<double>& v);
double sumCpp(std::vector<double>& v);
double meanCpp(std::vector<double>& v);
double maxCpp(std::vector<double>& v);
double minCpp(std::vector<double>& v);

bool lessCpp(double& x, double& y);
bool greaterCpp(double& x, double& y);
bool lessEqualCpp(double& x, double& y);
bool greaterEqualCpp(double& x, double& y);
bool equalCpp(double& x, double& y);

Rcpp::XPtr<funcPtr> putFunPtrInXPtr(std::string fstr);
Rcpp::XPtr<compPtr> putCompPtrInXPtr(std::string fstr);

#endif
