#ifndef CHECK_RETURN_H
#define CHECK_RETURN_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

extern "C" {
    SEXP CheckReturn(SEXP Rv, SEXP f1,
                     SEXP f2, SEXP Rtarget,
                     SEXP RKeepRes, SEXP stdFun);
}

bool CheckConstrnd(SEXP f1, SEXP f2, SEXP Rtarget);

#endif
