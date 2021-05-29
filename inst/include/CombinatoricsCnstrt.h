#ifndef COMBINATORICS_CNSTRT_H
#define COMBINATORICS_CNSTRT_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP CombinatoricsCnstrt(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                             SEXP Rlow, SEXP Rhigh, SEXP RmainFun,
                             SEXP RcompFun, SEXP Rtarget, SEXP RIsComb,
                             SEXP RKeepRes, SEXP Rparallel, SEXP RnThreads,
                             SEXP RmaxThreads, SEXP Rtolerance);
}

#endif
