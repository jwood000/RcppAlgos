#pragma once

#include "CppConvert.h"

SEXP GetCombPerms(SEXP Rv, const std::vector<double> &vNum,
                  const std::vector<int> &vInt, int n, int m, int phaseOne,
                  bool generalRet, bool IsComb, bool Parallel, bool IsRep,
                  bool IsMult, bool IsGmp, const std::vector<int> &freqs,
                  std::vector<int> &z, const std::vector<int> &myReps,
                  double lower, mpz_class &lowerMpz, int nRows,
                  int nThreads, VecType myType);
