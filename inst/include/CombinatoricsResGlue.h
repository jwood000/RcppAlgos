#ifndef COMBINATORICS_RES_GLUE_H
#define COMBINATORICS_RES_GLUE_H

#include "Permutations/PermResMain.h"
#include "Combinations/CombResMain.h"

template <typename T>
void ResultsMain(T* mat, const std::vector<T> &v, const funcPtr<T> myFun,
                 int n, int m, bool IsComb, bool Parallel, bool IsRep,
                 bool IsMult, bool IsGmp, const std::vector<int> &freqs,
                 std::vector<int> &z, const std::vector<int> &myReps,
                 double lower, mpz_t lowerMpz, int nRows, int nThreads) {
    
    if (IsComb) {
        CombResMain(mat, v, myFun, n, m, Parallel, IsRep,
                    IsMult, IsGmp, freqs, z, myReps, lower,
                    lowerMpz, nRows, nThreads);
    } else {
        PermResMain(mat, v, myFun, n, m, Parallel, IsRep,
                    IsMult, IsGmp, freqs, z, myReps, lower,
                    lowerMpz, nRows, nThreads);
    }
}

#endif
