#ifndef MAIN_PERMUTATIONS_RESULT_H
#define MAIN_PERMUTATIONS_RESULT_H

#include "Constraints/UserConstraintFuns.h"
#include <gmp.h>

template <typename T>
void PermuteResMain(T* mat, const std::vector<T> &v, const funcPtr<T> myFun,
                    int n, int m, bool Parallel, bool IsRep, bool IsMult,
                    bool IsGmp, const std::vector<int> &freqs,
                    std::vector<int> &z, const std::vector<int> &myReps,
                    double lower, mpz_t lowerMpz, int nRows, int nThreads);

#endif