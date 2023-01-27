#pragma once

#include "RMatrix.h"
#include <vector>

template <typename T>
void PartsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                 std::size_t width, int lastElem,
                 int lastCol, std::size_t nRows);

template <typename T>
void PartsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                 std::vector<int> &z, int strt, std::size_t width,
                 int lastElem, int lastCol, std::size_t nRows);

template <typename T>
void PartsGenPermRep(T*, const std::vector<T> &v, std::vector<int> &z,
                     std::size_t width, int lastElem,
                     int lastCol, std::size_t nRows);

template <typename T>
void PartsGenRep(std::vector<T> &partsVec, const std::vector<T> &v,
                 std::vector<int> &z, std::size_t width,
                 std::size_t nRows, bool IsComb);

void PartsRep(int* mat, std::vector<int> &z, std::size_t width,
              int lastElem, int lastCol, std::size_t nRows);

void PartsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
              int strt, std::size_t width, int lastElem,
              int lastCol, std::size_t nRows);

void PartsPermRep(int* mat, std::vector<int> &z, std::size_t width,
                  int lastElem, int lastCol, std::size_t nRows);
