#ifndef PERMUTE_RES_GLUE_H
#define PERMUTE_RES_GLUE_H

#include "Constraints/UserConstraintFuns.h"
#include "RMatrix.h"
#include <vector>

template <typename T>
void PermuteResStd(T* mat, const std::vector<T> &v,
                   std::vector<int> &z, int n, int m,
                   int nRows, bool IsMult, bool IsRep,
                   const std::vector<int> &freqs,
                   const funcPtr<T> myFun);

template <typename T>
void PermuteResPar(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                   std::vector<int> &z, int n, int m, int strt, int nRows,
                   const std::vector<int> &freqs, const funcPtr<T> myFun,
                   bool IsMult, bool IsRep);

#endif
