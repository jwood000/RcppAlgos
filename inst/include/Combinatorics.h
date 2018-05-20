#ifndef RcppAlgos_Combinatorics_h
#define RcppAlgos_Combinatorics_h

#include <Rcpp.h>

SEXP PermutationsRcpp(int n, int m, bool repetition, Rcpp::CharacterVector vStr,
                      int nRows, std::vector<int> vInt, std::vector<double> vNum,
                      bool isMult, bool isFac, bool keepRes, std::vector<int> z23,
                      bool isChar, SEXP Rv, bool isInt, std::vector<int> myReps,
                      SEXP f1, SEXP f2, bool nonTrivial);


SEXP CombinationsRcpp(int n, int m, bool repetition, Rcpp::CharacterVector vStr,
                      int nRows, std::vector<int> vInt, std::vector<double> vNum,
                      bool isMult, bool isFac, bool keepRes, std::vector<int> z,
                      bool isChar, SEXP Rv, bool isInt, std::vector<int> myReps,
                      SEXP f1, SEXP f2);

#endif
