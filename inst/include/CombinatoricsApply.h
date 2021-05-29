#ifndef COMBINATORICS_APPLY_H
#define COMBINATORICS_APPLY_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP CombinatoricsApply(SEXP Rv, SEXP Rm, SEXP RisRep,
                            SEXP RFreqs, SEXP Rlow, SEXP Rhigh,
                            SEXP stdFun, SEXP myEnv,
                            SEXP RFunVal, SEXP RIsComb);
}

#endif
