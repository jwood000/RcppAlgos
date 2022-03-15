#ifndef GET_PREV_COMB_PERM_H
#define GET_PREV_COMB_PERM_H

#include "ClassUtils/PrevCombinatorics.h"
#include "CleanConvert.h"

SEXP GetPrevCombPerms(SEXP Rv, const std::vector<double> &vNum,
                      const std::vector<int> &vInt,
                      const std::vector<int> &myReps,
                      const std::vector<int> &freqs, std::vector<int> &z,
                      prevIterPtr prevIter, int n, int m, bool IsComb,
                      bool IsMult, int nRows, VecType myType);

#endif
