#ifndef COMBINATORICS_MAIN_H
#define COMBINATORICS_MAIN_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP CombinatoricsStndrd(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                             SEXP Rlow, SEXP Rhigh, SEXP Rparallel,
                             SEXP RNumThreads, SEXP maxThreads, SEXP RIsComb);
}

#endif
