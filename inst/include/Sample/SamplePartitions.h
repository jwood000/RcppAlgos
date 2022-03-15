#ifndef SAMPLE_PARTITIONS_H
#define SAMPLE_PARTITIONS_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

#include "Partitions/NthPartition.h"

extern "C" {
    SEXP SamplePartitions(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                          SEXP RindexVec, SEXP RmySeed, SEXP RNumSamp,
                          SEXP baseSample, SEXP Rparallel, SEXP RNumThreads,
                          SEXP RmaxThreads, SEXP RNamed, SEXP RcompFun,
                          SEXP Rtarget, SEXP Rtolerance, SEXP myEnv);
}

template <typename T>
void ThreadSafeSample(T* mat, SEXP res, const std::vector<T> &v,
                      const std::vector<double> &mySample,
                      mpz_t *const myBigSamp, const std::vector<int> &myReps,
                      nthPartsPtr nthPartFun, int m, int sampSize,
                      int nThreads, bool Parallel, bool IsNamed,
                      int tar, int strtLen, int cap, bool IsGmp);

#endif
