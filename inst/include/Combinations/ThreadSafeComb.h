#pragma once

#include <vector>
#include "CppConvert/GmpxxCopy.h"

template <typename T>
void ThreadSafeCombinations(T* mat, const std::vector<T> &v, int n, int m,
                            bool Parallel, bool IsRep, bool IsMult, bool IsGmp,
                            const std::vector<int> &freqs, std::vector<int> &z,
                            const std::vector<int> &myReps, double lower,
                            mpz_class &lowerMpz, int nRows, int nThreads);
