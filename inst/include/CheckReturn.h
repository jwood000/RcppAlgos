#ifndef CHECK_RETURN_H
#define CHECK_RETURN_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP CheckReturn(SEXP Rv, SEXP RCnstrntFun, SEXP RCompFun,
                     SEXP Rtarget, SEXP RKeepRes, SEXP stdFun);
    SEXP CheckConstrndCpp(SEXP RCnstrntFun, SEXP RCompFun, SEXP Rtarget);
}

bool CheckConstrnd(SEXP RCnstrntFun, SEXP RCompFun, SEXP Rtarget);

#endif
