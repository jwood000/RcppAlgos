#ifndef PERMUTATION_RESULTS_H
#define PERMUTATION_RESULTS_H

#include "UserConstraintFuns.h"

template <typename typeMatrix, typename typeVector>
void PermuteGenResNoRep(typeMatrix &matRcpp, const std::vector<typeVector> &v, std::vector<int> z,
                        int n, int m, int strt, int nRows, const std::vector<int> &freqs,
                        funcPtr<typeVector> myFun);

template <typename typeMatrix, typename typeVector>
void PermuteGenResRep(typeMatrix &matRcpp, const std::vector<typeVector> &v, std::vector<int> z,
                      int n, int m, int strt, int nRows, const std::vector<int> &freqs,
                      funcPtr<typeVector> myFun);

template <typename typeMatrix, typename typeVector>
void MultisetPermRes(typeMatrix &matRcpp, const std::vector<typeVector> &v, std::vector<int> z,
                     int n, int m, int strt, int nRows, const std::vector<int> &freqs,
                     funcPtr<typeVector> myFun);

#endif
