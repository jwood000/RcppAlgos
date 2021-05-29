#ifndef COMBINATORICS_GLUE_H
#define COMBINATORICS_GLUE_H

#include "Permutations/ThreadSafePerm.h"
#include "Combinations/ThreadSafeComb.h"
#include "Permutations/PermuteManager.h"
#include "Combinations/ComboManager.h"

void CharacterGlue(SEXP mat, SEXP v, bool IsComb,
                   std::vector<int> &z, int n, int m, int nRows,
                   const std::vector<int> &freqs, bool IsMult, bool IsRep) {

    if (IsComb) {
        ComboCharacter(mat, v, z, n, m, nRows, freqs, IsMult, IsRep);
    } else {
        PermuteCharacter(mat, v, z, n, m, nRows, freqs, IsMult, IsRep);
    }
}

template <typename T>
void ManagerGlue(T* mat, const std::vector<T> &v, std::vector<int> &z,
                 int n, int m, int nRows, bool IsComb, int phaseOne,
                 bool generalRet, const std::vector<int> &freqs,
                 bool IsMult, bool IsRep) {

    if (IsComb) {
        ComboManager(mat, v, z, n, m, nRows, freqs, IsMult, IsRep);
    } else {
        PermuteManager(mat, v, z, n, m, nRows, phaseOne,
                       generalRet, IsMult, IsRep, freqs);
    }
}

template <typename T>
void ParallelGlue(T* mat, const std::vector<T> &v, int n, int m, int phaseOne,
                  bool generalRet, bool IsComb, bool Parallel, bool IsRep,
                  bool IsMult, bool IsGmp, const std::vector<int> &freqs,
                  std::vector<int> &z, const std::vector<int> &myReps,
                  double lower, mpz_t lowerMpz, int nRows, int nThreads) {

    if (IsComb) {
        ThreadSafeCombinations(mat, v, n, m, Parallel, IsRep,
                               IsMult, IsGmp, freqs, z, myReps,
                               lower, lowerMpz, nRows, nThreads);
    } else {
        ThreadSafePermutations(mat, v, n, m, phaseOne, generalRet, Parallel,
                               IsRep, IsMult, IsGmp, freqs, z, myReps, lower,
                               lowerMpz, nRows, nThreads);
    }
}

#endif
