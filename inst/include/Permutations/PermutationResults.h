#pragma once

#include "Permutations/NextPermSectionRep.h"
#include "Constraints/UserConstraintFuns.h"
#include "Permutations/NextPermutation.h"
#include <memory>
#include "RMatrix.h"

template <typename T>
void PermuteResDistinct(T* mat, const std::vector<T> &v,
                        std::vector<int> &z, std::size_t n, std::size_t m,
                        std::size_t nRows, const funcPtr<T> myFun);

template <typename T>
void PermuteResDistinct(RcppParallel::RMatrix<T> &mat,
                        const std::vector<T> &v, std::vector<int> &z,
                        std::size_t n, std::size_t m, std::size_t strt,
                        std::size_t nRows, const funcPtr<T> myFun);

template <typename T>
void PermuteResRep(T* mat, const std::vector<T> &v,
                   std::vector<int> &z, std::size_t n, std::size_t m,
                   std::size_t nRows, const funcPtr<T> myFun);

template <typename T>
void PermuteResRep(RcppParallel::RMatrix<T> &mat,
                   const std::vector<T> &v,
                   std::vector<int> &z, std::size_t n,
                   std::size_t m, std::size_t strt, std::size_t nRows,
                   const funcPtr<T> myFun);

template <typename T>
void MultisetPermRes(T* mat, const std::vector<T> &v,
                     std::vector<int> &z, std::size_t n, std::size_t m,
                     std::size_t nRows, const std::vector<int> &freqs,
                     const funcPtr<T> myFun);

template <typename T>
void MultisetPermRes(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     std::vector<int> &z, std::size_t n, std::size_t m,
                     std::size_t strt, std::size_t nRows,
                     const std::vector<int> &freqs, const funcPtr<T> myFun);
