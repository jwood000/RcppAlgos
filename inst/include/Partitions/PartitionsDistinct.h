#pragma once

#include "RMatrix.h"
#include <vector>

template <typename T>
int PartsGenDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, std::size_t width,
                      int lastElem, int lastCol, std::size_t nRows);

template <typename T>
int PartsGenDistinct(RcppParallel::RMatrix<T> &mat,
                     const std::vector<T> &v, std::vector<int> &z,
                     int strt, std::size_t width, int lastElem,
                     int lastCol, std::size_t nRows);

template <typename T>
int PartsGenDistinct(std::vector<T> &partsVec, const std::vector<T> &v,
                     std::vector<int> &z, std::size_t width,
                     std::size_t nRows, bool IsComb);

template <typename T>
int PartsGenPermDistinct(T* mat, const std::vector<T> &v,
                         std::vector<int> &z, std::size_t width,
                         int lastElem, int lastCol, std::size_t nRows);

template <typename T>
int PartsGenPermZeroDistinct(T* mat, const std::vector<T> &v,
                             std::vector<int> &z, std::size_t width,
                             int lastElem, int lastCol, std::size_t nRows);

int PartsDistinct(int* mat, std::vector<int> &z, std::size_t width,
                  int lastElem, int lastCol, std::size_t nRows);

int PartsDistinct(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                  int strt, std::size_t width, int lastElem,
                  int lastCol, std::size_t nRows);

int PartsPermDistinct(int* mat, std::vector<int> &z, std::size_t width,
                      int lastElem, int lastCol, std::size_t nRows);

int PartsPermZeroDistinct(int* mat, std::vector<int> &z, std::size_t width,
                          int lastElem, int lastCol, std::size_t nRows);
