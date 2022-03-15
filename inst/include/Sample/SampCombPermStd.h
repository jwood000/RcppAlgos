#ifndef SAMP_COMB_PERM_STD_H
#define SAMP_COMB_PERM_STD_H

#include "CleanConvert.h"
#include "NthResult.h"

SEXP SampCombPermMain(SEXP Rv, const std::vector<int> &vInt,
                      const std::vector<double> &vNum,
                      const std::vector<double> &mySample,
                      mpz_t *const myBigSamp,
                      const std::vector<int> &myReps,
                      nthResultPtr nthResFun, VecType myType, int n,
                      int m, int sampSize, int nThreads, bool IsNamed,
                      bool IsGmp, bool Parallel);

#endif
