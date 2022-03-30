#ifndef COMBO_MANAGER_H
#define COMBO_MANAGER_H

#include "cpp11/R.hpp"
#include <vector>
#include "RMatrix.h"

template <typename T>
void ComboManager(T* mat, const std::vector<T> &v,
                  std::vector<int> &z, int n, int m, int nRows,
                  const std::vector<int> &freqs, bool IsMult, bool IsRep);

template <typename T>
void ComboParallel(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                   std::vector<int> &z, int n, int m, int strt, int nRows,
                   const std::vector<int> &freqs, bool IsMult, bool IsRep);

void ComboCharacter(SEXP mat, SEXP v, std::vector<int> &z, int n,
                    int m, int nRows, const std::vector<int> &freqs,
                    bool IsMult, bool IsRep);

#endif
