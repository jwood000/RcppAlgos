#ifndef GET_COMB_PERM_H
#define GET_COMB_PERM_H

#include "CleanConvert.h"

SEXP GetCombPerms(SEXP Rv, const std::vector<double> &vNum,
                  const std::vector<int> &vInt, int n, int m, int phaseOne,
                  bool generalRet, bool IsComb, bool Parallel, bool IsRep,
                  bool IsMult, bool IsGmp, const std::vector<int> &freqs,
                  std::vector<int> &z, const std::vector<int> &myReps,
                  double lower, mpz_t lowerMpz, int nRows,
                  int nThreads, VecType myType);

#endif
