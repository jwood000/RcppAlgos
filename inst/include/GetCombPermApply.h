#pragma once

#include "CppConvert.h"

SEXP GetCombPermApply(SEXP Rv, const std::vector<double> &vNum,
                      const std::vector<int> &vInt, int n, int m,
                      bool IsComb, bool IsRep, bool IsMult,
                      const std::vector<int> &freqs, std::vector<int> &z,
                      const std::vector<int> &myReps, VecType myType,
                      int nRows, SEXP stdFun, SEXP myEnv, SEXP RFunVal);
