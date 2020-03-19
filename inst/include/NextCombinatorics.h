#ifndef NEXT_COMBINATORICS_H
#define NEXT_COMBINATORICS_H

#include <Rcpp.h>

using nextIterPtr = void (*const)(const std::vector<int> &freqs,
                          std::vector<int> &z,int n1, int m1);

Rcpp::XPtr<nextIterPtr> putNextIterPtrInXPtr(bool IsComb, bool IsMult, bool IsRep, bool IsFull);

#endif
