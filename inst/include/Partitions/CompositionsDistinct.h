#pragma once

#include "RMatrix.h"
#include <vector>

template <typename T>
void CompsGenDistinct(
    T* mat, const std::vector<T> &v, std::vector<int> &z,
    std::vector<int> &complement, int i1, int i2, int myMax,
    int tar, std::size_t width, std::size_t nRows
);

template <typename T>
void CompsGenDistinct(
    RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
    std::vector<int> &z, std::vector<int> &complement, int i1, int i2,
    int myMax, int tar, std::size_t strt, std::size_t width, std::size_t nRows
);

void CompsDistinct(int* mat, std::vector<int> &z,
                   std::size_t width, std::size_t nRows);

void CompsDistinct(
    RcppParallel::RMatrix<int> &mat,
    std::vector<int> &z, std::vector<int> &complement, int i1, int i2,
    int myMax, int tar, std::size_t strt, std::size_t width, std::size_t nRows
);
