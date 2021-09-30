#ifndef SAMPLE_APPLY_H
#define SAMPLE_APPLY_H

#include "CleanConvert.h"
#include "NthResult.h"

SEXP SampleCombPermApply(SEXP Rv, const std::vector<int> &vInt,
                         const std::vector<double> &vNum,
                         const std::vector<double> &mySample,
                         mpz_t *const myBigSamp,
                         const std::vector<int> &myReps, SEXP stdFun,
                         SEXP rho, SEXP RFunVal, nthResultPtr nthResFun,
                         VecType myType, int n, int m, int sampSize,
                         bool IsNamed, bool IsGmp);

#endif
