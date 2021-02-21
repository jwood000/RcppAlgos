#ifndef COMBINATION_RESULTS_H
#define COMBINATION_RESULTS_H

#include "Constraints/UserConstraintFuns.h"
#include "Combinations/NextComboSection.h"
#include "RMatrix.h"

template <typename T>
void ComboGenResDistinct(T* mat, const std::vector<T> &v,
                         std::vector<int> &z, int n, int m, int strt,
                         int nRows, const funcPtr<T> myFun);
template <typename T>
void ComboGenResDistinct(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                         std::vector<int> &z, int n, int m, int strt,
                         int nRows, const funcPtr<T> myFun);
template <typename T>
void ComboGenResRep(T* mat, const std::vector<T> &v, 
                    std::vector<int> &z, int n, int m, int strt,
                    int nRows, const funcPtr<T> myFun);
template <typename T>
void ComboGenResRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v, 
                    std::vector<int> &z, int n, int m, int strt,
                    int nRows, const funcPtr<T> myFun);
template <typename T>
void MultisetComboResult(T* mat, const std::vector<T> &v,
                         std::vector<int> &z, int n, int m, int strt, int nRows,
                         const std::vector<int> &freqs, const funcPtr<T> myFun);
template <typename T>
void MultisetComboResult(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                         std::vector<int> &z, int n, int m, int strt, int nRows,
                         const std::vector<int> &freqs, const funcPtr<T> myFun);

#endif
