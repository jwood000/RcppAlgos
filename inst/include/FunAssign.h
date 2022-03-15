#ifndef FUN_ASSIGN_H
#define FUN_ASSIGN_H

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

void FunAssign(SEXP res, SEXP vectorPass, SEXP sexpFun,
               SEXP rho, int commonType, int commonLen,
               int count, int nRows, int retType);

void SetDims(SEXP RFunVal, SEXP res, int commonLen, int nRows);

#endif
