#ifndef COMPOSITIONS_REP_H
#define COMPOSITIONS_REP_H

#include "RMatrix.h"
#include <vector>

template <int one_or_zero, typename T>
void CompsGenRep(T* mat, const std::vector<T> &v,
                 std::vector<int> &z, int width, int nRows);

template <int one_or_zero, typename T>
void CompsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                 std::vector<int> &z, int strt, int width, int nRows);

template <int one_or_zero>
void CompsRep(int* mat, std::vector<int> &z,
              int width, int nRows);

template <int one_or_zero>
void CompsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
              int strt, int width, int nRows);

#endif
