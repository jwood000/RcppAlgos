#ifndef THREAD_SAFE_COMBINATIONS_H
#define THREAD_SAFE_COMBINATIONS_H

#include <vector>
#include <gmp.h>

template <typename T>
void ThreadSafeCombinations(T* mat, const std::vector<T> &v, int n, int m,
                            bool Parallel, bool IsRep, bool IsMult, bool IsGmp,
                            const std::vector<int> &freqs, std::vector<int> &z,
                            const std::vector<int> &myReps, double lower,
                            mpz_t lowerMpz, int nRows, int nThreads);

#endif