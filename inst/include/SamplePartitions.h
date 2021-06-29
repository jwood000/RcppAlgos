#ifndef SAMPLE_PARTITIONS_H
#define SAMPLE_PARTITIONS_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP SamplePartitions(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                          SEXP RindexVec, SEXP RmySeed, SEXP RNumSamp,
                          SEXP baseSample, SEXP Rparallel, SEXP RNumThreads,
                          SEXP RmaxThreads, SEXP RNamed, SEXP RcompFun,
                          SEXP Rtarget, SEXP Rtolerance, SEXP myEnv);
}

#endif
