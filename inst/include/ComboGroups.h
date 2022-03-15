#ifndef COMBO_GROUPS_H
#define COMBO_GROUPS_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP ComboGroupsCpp(SEXP Rv, SEXP RNumGroups, SEXP RRetType, SEXP Rlow,
                        SEXP Rhigh, SEXP Rparallel, SEXP RNumThreads,
                        SEXP RmaxThreads, SEXP RIsSample, SEXP RindexVec,
                        SEXP RmySeed, SEXP RNumSamp, SEXP baseSample,
                        SEXP RNamed, SEXP myEnv);
}

#endif
