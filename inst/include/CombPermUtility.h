#ifndef RcppAlgos_CombPermUtility_h
#define RcppAlgos_CombPermUtility_h

#include <Rcpp.h>

typedef double (*funcPtr)(std::vector<double>& x);
typedef bool (*compPtr)(double& x, double& y);

double NumPermsWithRep(std::vector<int> v);
double NumPermsNoRep(int n, int k);
double nChooseK(double n, double k);
double NumCombsWithRep(int n, int r);
double MultisetCombRowNum(int n, int r, std::vector<int> Reps);
double MultisetPermRowNum(int n, int r, std::vector<int> Reps, 
                          Rcpp::IntegerMatrix myCombs);

// Contraint functions
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
