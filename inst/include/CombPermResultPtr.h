#ifndef COMB_PERM_RESULT_PTR_H
#define COMB_PERM_RESULT_PTR_H

#include "CombinationResults.h"
#include "PermutationResults.h"

template <typename T, typename U>
using combPermResPtr = void (*const)(T &matRcpp, const std::vector<U> &v,
                             std::vector<int> z, int n, int m, int strt, int nRows,
                             const std::vector<int> &freqs, funcPtr<U> myFun);

template <typename T, typename U>
Rcpp::XPtr<combPermResPtr<T, U>> putCombResPtrInXPtr(bool IsComb, bool IsMult, bool IsRep) {
    if (IsComb) {
        if (IsMult)
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&MultisetComboResult)));
        else if (IsRep)
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&ComboGenResRep)));
        else
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&ComboGenResNoRep)));
    } else {
        if (IsMult)
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&MultisetPermRes)));
        else if (IsRep)
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&PermuteGenResRep)));
        else
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&PermuteGenResNoRep)));
    }
}

#endif
