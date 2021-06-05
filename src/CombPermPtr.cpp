#include "Combinations.h"
#include "Permutations.h"
#include "RMatrix.h"

template <typename T, typename U>
using combPermPtr = void (*const)(T &matRcpp, const U &v, std::vector<int> z,
                          int n, int m, int strt, int nRows, const std::vector<int> &freqs);

template <typename T, typename U>
Rcpp::XPtr<combPermPtr<T, U>> putCombPtrInXPtr(bool IsComb, bool IsMult, bool IsRep, bool IsGen) {
    
    if (IsComb) {
        if (IsMult)
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&MultisetCombination)));
        else if (IsRep)
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&CombinationsRep)));
        else
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&CombinationsNoRep)));
    } else {
        if (IsMult) {
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&MultisetPermutation)));
        } else if (IsGen) {
            if (IsRep)
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteGeneralRep)));
            else
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteGeneralNoRep)));
        } else {
            if (IsRep)
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteSerialRep)));
            else
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteSerialNoRep)));
        }
    }
}

template Rcpp::XPtr<combPermPtr<Rcpp::IntegerMatrix,
                                std::vector<int>>> putCombPtrInXPtr(bool, bool, bool, bool);

template Rcpp::XPtr<combPermPtr<Rcpp::NumericMatrix,
                                std::vector<double>>> putCombPtrInXPtr(bool, bool, bool, bool);

template Rcpp::XPtr<combPermPtr<RcppParallel::RMatrix<int>,
                                std::vector<int>>> putCombPtrInXPtr(bool, bool, bool, bool);

template Rcpp::XPtr<combPermPtr<RcppParallel::RMatrix<double>,
                                std::vector<double>>> putCombPtrInXPtr(bool, bool, bool, bool);

template Rcpp::XPtr<combPermPtr<Rcpp::IntegerMatrix,
                                Rcpp::IntegerVector>> putCombPtrInXPtr(bool, bool, bool, bool);

template Rcpp::XPtr<combPermPtr<Rcpp::NumericMatrix,
                                Rcpp::NumericVector>> putCombPtrInXPtr(bool, bool, bool, bool);

template Rcpp::XPtr<combPermPtr<Rcpp::LogicalMatrix,
                                Rcpp::LogicalVector>> putCombPtrInXPtr(bool, bool, bool, bool);

template Rcpp::XPtr<combPermPtr<Rcpp::CharacterMatrix,
                                Rcpp::CharacterVector>> putCombPtrInXPtr(bool, bool, bool, bool);

template Rcpp::XPtr<combPermPtr<Rcpp::ComplexMatrix,
                                Rcpp::ComplexVector>> putCombPtrInXPtr(bool, bool, bool, bool);

template Rcpp::XPtr<combPermPtr<Rcpp::RawMatrix,
                                Rcpp::RawVector>> putCombPtrInXPtr(bool, bool, bool, bool);

    
