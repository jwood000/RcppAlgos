#pragma once

#include "RMatrix.h"
#include <vector>

template <typename T>
int PartsGenMultiset(T* mat, const std::vector<T> &v,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     std::size_t width, int lastElem,
                     int lastCol, std::size_t nRows);

template <typename T>
int PartsGenMultiset(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     int strt, std::size_t width, int lastElem,
                     int lastCol, std::size_t nRows);

template <typename T>
int PartsGenPermMultiset(T* mat, const std::vector<T> &v,
                         const std::vector<int> &Reps, std::vector<int> &z,
                         std::size_t width, int lastElem,
                         int lastCol, std::size_t nRows);
