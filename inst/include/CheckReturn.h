#ifndef CHECK_RETURN_H
#define CHECK_RETURN_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP CheckReturn(SEXP Rv, SEXP RCnstrntFun, SEXP RCompFun,
                     SEXP Rtarget, SEXP RKeepRes, SEXP stdFun);
    
    SEXP CheckPartition(SEXP Rv, SEXP Rm, SEXP RmainFun, SEXP RcompFun,
                        SEXP Rlow, SEXP Rtarget, SEXP Rtolerance);
}

bool CheckConstrnd(SEXP f1, SEXP f2, SEXP Rtarget);

#endif
