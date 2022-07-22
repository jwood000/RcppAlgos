#ifndef COMPOSITIONS_REP_H
#define COMPOSITIONS_REP_H

#include "RMatrix.h"
#include <vector>

template <typename T>
void CompsGenRep(T* mat, const std::vector<T> &v,
                 std::vector<int> &z, int width, int nRows);

template <typename T>
void CompsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                 std::vector<int> &z, int strt, int width, int nRows);

void CompsRep(int* mat, std::vector<int> &z,
              int width, int nRows);

void CompsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
              int strt, int width, int nRows);

#endif
