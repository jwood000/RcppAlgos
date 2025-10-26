#pragma once

#include "RMatrix.h"
#include <vector>

template <typename T>
int PartsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                std::size_t width, int lastElem,
                int lastCol, std::size_t nRows);

template <typename T>
int PartsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                std::vector<int> &z, int strt, std::size_t width,
                int lastElem, int lastCol, std::size_t nRows);

template <typename T>
int PartsGenPermRep(T*, const std::vector<T> &v, std::vector<int> &z,
                    std::size_t width, int lastElem,
                    int lastCol, std::size_t nRows);

template <typename T>
int PartsGenRep(std::vector<T> &partsVec, const std::vector<T> &v,
                std::vector<int> &z, std::size_t width,
                std::size_t nRows, bool IsComb);

int PartsRep(int* mat, std::vector<int> &z, std::size_t width,
             int lastElem, int lastCol, std::size_t nRows);

int PartsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
             int strt, std::size_t width, int lastElem,
             int lastCol, std::size_t nRows);

int PartsPermRep(int* mat, std::vector<int> &z, std::size_t width,
                 int lastElem, int lastCol, std::size_t nRows);
