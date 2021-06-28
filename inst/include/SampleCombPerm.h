#ifndef SAMPLE_COMB_PERM_H
#define SAMPLE_COMB_PERM_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP SampleCombPerm(
            SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP RindexVec,
            SEXP RIsComb, SEXP RmySeed, SEXP RNumSamp, SEXP baseSample,
            SEXP stdFun, SEXP myEnv, SEXP Rparallel, SEXP RNumThreads,
            SEXP RmaxThreads, SEXP RNamed, SEXP RFunVal
    );
}

#endif
