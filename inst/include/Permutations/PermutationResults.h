#ifndef PERMUTATION_RESULTS_H
#define PERMUTATION_RESULTS_H

#include "Permutations/NextPermSectionRep.h"
#include "Constraints/UserConstraintFuns.h"
#include "Permutations/NextPermutation.h"
#include "Cpp14MakeUnique.h"
#include "RMatrix.h"

template <typename T>
void PermuteResDistinct(T* mat, const std::vector<T> &v,
                        std::vector<int> &z, int n, int m,
                        int nRows, const funcPtr<T> myFun);

template <typename T>
void PermuteResDistinct(RcppParallel::RMatrix<T> &mat,
                        const std::vector<T> &v, std::vector<int> &z,
                        int n, int m, int strt, int nRows,
                        const funcPtr<T> myFun);

template <typename T>
void PermuteResRep(T* mat, const std::vector<T> &v,
                   std::vector<int> &z, int n, int m,
                   int nRows, const funcPtr<T> myFun);

template <typename T>
void PermuteResRep(RcppParallel::RMatrix<T> &mat,
                   const std::vector<T> &v,
                   std::vector<int> &z, int n, int m,
                   int strt, int nRows, const funcPtr<T> myFun);

template <typename T>
void MultisetPermRes(T* mat, const std::vector<T> &v,
                     std::vector<int> &z, int n, int m, int nRows,
                     const std::vector<int> &freqs, const funcPtr<T> myFun);

template <typename T>
void MultisetPermRes(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     std::vector<int> &z, int n, int m, int strt, int nRows,
                     const std::vector<int> &freqs, const funcPtr<T> myFun);

#endif
