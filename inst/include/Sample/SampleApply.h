#pragma once

#include "CppConvert.h"
#include "NthResult.h"

SEXP SampleCombPermApply(SEXP Rv, const std::vector<int> &vInt,
                         const std::vector<double> &vNum,
                         const std::vector<double> &mySample,
                         const std::vector<mpz_class> &myBigSamp,
                         const std::vector<int> &myReps, SEXP stdFun,
                         SEXP rho, SEXP RFunVal, nthResultPtr nthResFun,
                         VecType myType, int n, int m, int sampSize,
                         bool IsNamed, bool IsGmp);
