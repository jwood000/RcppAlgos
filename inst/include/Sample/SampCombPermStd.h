#pragma once

#include "CppConvert.h"
#include "NthResult.h"

SEXP SampCombPermMain(SEXP Rv, const std::vector<int> &vInt,
                      const std::vector<double> &vNum,
                      const std::vector<double> &mySample,
                      const std::vector<mpz_class> &myBigSamp,
                      const std::vector<int> &myReps,
                      nthResultPtr nthResFun, VecType myType, int n,
                      int m, int sampSize, int nThreads, bool IsNamed,
                      bool IsGmp, bool Parallel);
