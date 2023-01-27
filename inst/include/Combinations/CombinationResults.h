#pragma once

#include "Constraints/UserConstraintFuns.h"
#include "Combinations/NextComboSection.h"
#include "RMatrix.h"

template <typename T>
void ComboResDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, std::size_t n, std::size_t m,
                      std::size_t nRows, const funcPtr<T> myFun);
template <typename T>
void ComboResDistinct(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                      std::vector<int> &z, std::size_t n, std::size_t m,
                      std::size_t strt, std::size_t nRows,
                      const funcPtr<T> myFun);
template <typename T>
void ComboResRep(T* mat, const std::vector<T> &v,
                 std::vector<int> &z, std::size_t n, std::size_t m,
                 std::size_t nRows, const funcPtr<T> myFun);

template <typename T>
void ComboResRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                 std::vector<int> &z, std::size_t n, std::size_t m,
                 std::size_t strt, std::size_t nRows, const funcPtr<T> myFun);

template <typename T>
void MultisetComboResult(T* mat, const std::vector<T> &v,
                         std::vector<int> &z, std::size_t n, std::size_t m,
                         std::size_t nRows, const std::vector<int> &freqs,
                         const funcPtr<T> myFun);
template <typename T>
void MultisetComboResult(RcppParallel::RMatrix<T> &mat,
                         const std::vector<T> &v, std::vector<int> &z,
                         std::size_t n, std::size_t m, std::size_t strt,
                         std::size_t nRows, const std::vector<int> &freqs,
                         const funcPtr<T> myFun);
