#ifndef PREV_COMBINATORICS_H
#define PREV_COMBINATORICS_H

#include <Rcpp.h>

using prevIterPtr = void (*const)(const std::vector<int> &freqs,
                          std::vector<int> &z,int n1, int m1);

Rcpp::XPtr<prevIterPtr> putPrevIterPtrInXPtr(bool IsComb, bool IsMult, bool IsRep, bool IsFull);

#endif
