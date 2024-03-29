#pragma once

#include "Constraints/UserConstraintFuns.h"
#include <gmpxx.h>

template <typename T>
void PermuteResMain(T* mat, const std::vector<T> &v, const funcPtr<T> myFun,
                    int n, int m, bool Parallel, bool IsRep, bool IsMult,
                    bool IsGmp, const std::vector<int> &freqs,
                    std::vector<int> &z, const std::vector<int> &myReps,
                    double lower, mpz_class &lowerMpz,
                    int nRows, int nThreads);
