#pragma once

#include "ClassUtils/PrevCombinatorics.h"
#include "CppConvert.h"

SEXP GetPrevCombPermApply(SEXP Rv, const std::vector<double> &vNum,
                          const std::vector<int> &vInt,
                          const std::vector<int> &myReps,
                          const std::vector<int> &freqs, std::vector<int> &z,
                          prevIterPtr prevIter, int n, int m, bool IsComb,
                          bool IsMult, int nRows, VecType myType,
                          SEXP stdFun, SEXP myEnv, SEXP RFunVal);
