#ifndef FUN_ASSIGN_H
#define FUN_ASSIGN_H

#include "cpp11/R.hpp"

void FunAssign(SEXP res, SEXP vectorPass, SEXP sexpFun,
               SEXP rho, int commonType, int commonLen,
               int count, int nRows, int retType);

void SetDims(SEXP RFunVal, SEXP res, int commonLen, int nRows);

#endif
