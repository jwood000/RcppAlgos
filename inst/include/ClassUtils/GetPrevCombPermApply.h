#ifndef GET_PREV_COMB_PERM_APPLY_H
#define GET_PREV_COMB_PERM_APPLY_H

#include "ClassUtils/PrevCombinatorics.h"
#include "CleanConvert.h"

SEXP GetPrevCombPermApply(SEXP Rv, const std::vector<double> &vNum,
                          const std::vector<int> &vInt,
                          const std::vector<int> &myReps,
                          const std::vector<int> &freqs, std::vector<int> &z,
                          prevIterPtr prevIter, int n, int m, bool IsComb,
                          bool IsMult, int nRows, VecType myType,
                          SEXP stdFun, SEXP myEnv, SEXP RFunVal);

#endif
