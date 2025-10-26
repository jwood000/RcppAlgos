#pragma once

#include "RMatrix.h"
#include <vector>

template <typename T>
int CompsGenDistinct(T* mat, const std::vector<T> &v, std::vector<int> &z,
                     std::size_t width, std::size_t nRows);

template <typename T>
int CompsGenDistinct(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     std::vector<int> &z, std::size_t strt,
                     std::size_t width, std::size_t nRows);

int CompsDistinct(int* mat, std::vector<int> &z,
                   std::size_t width, std::size_t nRows);

int CompsDistinct(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                  std::size_t strt, std::size_t width, std::size_t nRows);
