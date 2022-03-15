#ifndef COMBINATORICS_COUNT_H
#define COMBINATORICS_COUNT_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP CombinatoricsCount(SEXP Rv, SEXP Rm, SEXP RisRep,
                            SEXP RFreqs, SEXP RIsComb);
    SEXP PartitionsCount(SEXP Rtarget, SEXP Rv, SEXP Rm,
                         SEXP RisRep, SEXP RFreqs, SEXP RcompFun,
                         SEXP Rlow, SEXP Rtolerance,
                         SEXP RPartDesign, SEXP Rshow);
    SEXP ComboGroupsCountCpp(SEXP Rv, SEXP RNumGroups);
}

#endif
