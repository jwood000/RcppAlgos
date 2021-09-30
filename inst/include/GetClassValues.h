#ifndef GET_CLASS_VALUES_H
#define GET_CLASS_VALUES_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP GetClassVals(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                      SEXP RIsComb, SEXP stdFun, SEXP RThreads,
                      SEXP RmaxThreads, SEXP RIsCnstrd);
}

#endif
