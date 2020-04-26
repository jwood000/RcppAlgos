#ifndef COMB_PERM_PTR_H
#define COMB_PERM_PTR_H


template <typename T, typename U>
using combPermPtr = void (*const)(T &matRcpp, const U &v, std::vector<int> z,
                          int n, int m, int strt, int nRows, const std::vector<int> &freqs);

template <typename T, typename U>
Rcpp::XPtr<combPermPtr<T, U>> putCombPtrInXPtr(bool IsComb, bool IsMult, bool IsRep, bool IsGen);

#endif
