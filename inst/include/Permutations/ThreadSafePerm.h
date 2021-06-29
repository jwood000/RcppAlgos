#ifndef THREAD_SAFE_PERMUTATIONS_H
#define THREAD_SAFE_PERMUTATIONS_H

#include <gmp.h>

template <typename T>
void ThreadSafePermutations(T* mat, const std::vector<T> &v, int n, int m,
                            int phaseOne, bool generalRet, bool Parallel,
                            bool IsRep, bool IsMult, bool IsGmp,
                            const std::vector<int> &freqs,
                            std::vector<int> &z,
                            const std::vector<int> &myReps, double lower,
                            mpz_t lowerMpz, int nRows, int nThreads);

#endif