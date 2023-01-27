#pragma once

#include "RMatrix.h"
#include <vector>

template <typename T>
void PartsGenDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, std::size_t width,
                      int lastElem, int lastCol, std::size_t nRows);

template <typename T>
void PartsGenDistinct(RcppParallel::RMatrix<T> &mat,
                      const std::vector<T> &v, std::vector<int> &z,
                      int strt, std::size_t width, int lastElem,
                      int lastCol, std::size_t nRows);

template <typename T>
void PartsGenDistinct(std::vector<T> &partsVec, const std::vector<T> &v,
                      std::vector<int> &z, std::size_t width,
                      std::size_t nRows, bool IsComb);

template <typename T>
void PartsGenPermDistinct(T* mat, const std::vector<T> &v,
                          std::vector<int> &z, std::size_t width,
                          int lastElem, int lastCol, std::size_t nRows);

template <typename T>
void PartsGenPermZeroDistinct(T* mat, const std::vector<T> &v,
                              std::vector<int> &z, std::size_t width,
                              int lastElem, int lastCol, std::size_t nRows);

void PartsDistinct(int* mat, std::vector<int> &z, std::size_t width,
                   int lastElem, int lastCol, std::size_t nRows);

void PartsDistinct(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                   int strt, std::size_t width, int lastElem,
                   int lastCol, std::size_t nRows);

void PartsPermDistinct(int* mat, std::vector<int> &z, std::size_t width,
                       int lastElem, int lastCol, std::size_t nRows);

void PartsPermZeroDistinct(int* mat, std::vector<int> &z, std::size_t width,
                           int lastElem, int lastCol, std::size_t nRows);
