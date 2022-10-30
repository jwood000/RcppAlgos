#pragma once

#include "Permutations/PermuteResMain.h"
#include "Combinations/ComboResMain.h"

template <typename T>
void ResultsMain(T* mat, const std::vector<T> &v, const funcPtr<T> myFun,
                 int n, int m, bool IsComb, bool Parallel, bool IsRep,
                 bool IsMult, bool IsGmp, const std::vector<int> &freqs,
                 std::vector<int> &z, const std::vector<int> &myReps,
                 double lower, mpz_class &lowerMpz, int nRows, int nThreads) {

    if (IsComb) {
        ComboResMain(mat, v, myFun, n, m, Parallel, IsRep,
                     IsMult, IsGmp, freqs, z, myReps, lower,
                     lowerMpz, nRows, nThreads);
    } else {
        PermuteResMain(mat, v, myFun, n, m, Parallel, IsRep,
                       IsMult, IsGmp, freqs, z, myReps, lower,
                       lowerMpz, nRows, nThreads);
    }
}
