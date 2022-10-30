#pragma once

#include "cpp11/R.hpp"
#include "Partitions/NthPartition.h"

template <typename T>
void ThreadSafeSample(T* mat, SEXP res, const std::vector<T> &v,
                      const std::vector<double> &mySample,
                      const std::vector<mpz_class> &myBigSamp,
                      const std::vector<int> &myReps,
                      nthPartsPtr nthPartFun, int m, int sampSize,
                      int nThreads, bool Parallel, bool IsNamed,
                      int tar, int strtLen, int cap, bool IsGmp);
