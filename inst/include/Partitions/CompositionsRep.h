#pragma once

#include "RMatrix.h"
#include <vector>

template <int one_or_zero, typename T>
int CompsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                std::size_t width, std::size_t nRows);

template <int one_or_zero, typename T>
int CompsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                std::vector<int> &z, std::size_t strt,
                std::size_t width, std::size_t nRows);

template <int one_or_zero>
int CompsRep(int* mat, std::vector<int> &z,
             std::size_t width, std::size_t nRows);

template <int one_or_zero>
int CompsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
             std::size_t strt, std::size_t width, std::size_t nRows);
