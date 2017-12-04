#ifndef RcppAlgos_CombPermUtility_h
#define RcppAlgos_CombPermUtility_h

double NumPermsWithRep(std::vector<double> v);
double NumPermsNoRep(int n, int k);
double nChooseK(double n, double k);
double GetRowNum(int n, int r);

// Contraint functions
double prodCpp(std::vector<double>& v);
double sumCpp(std::vector<double>& v);
double meanCpp(std::vector<double>& v);
double maxCpp(std::vector<double>& v);
double minCpp(std::vector<double>& v);

#endif