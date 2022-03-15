#ifndef PARTITIONS_DISTINCT_H
#define PARTITIONS_DISTINCT_H

#include "RMatrix.h"
#include <vector>

template <typename T>
void PartsGenDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, int width,
                      int lastElem, int lastCol, int nRows);

template <typename T>
void PartsGenDistinct(RcppParallel::RMatrix<T> &mat,
                      const std::vector<T> &v, std::vector<int> &z,
                      int strt, int width, int lastElem,
                      int lastCol, int nRows);

template <typename T>
void PartsGenDistinct(std::vector<T> &partsVec, const std::vector<T> &v,
                      std::vector<int> &z, int width,
                      int nRows, bool IsComb);

template <typename T>
void PartsGenPermDistinct(T* mat, const std::vector<T> &v,
                          std::vector<int> &z, int width,
                          int lastElem, int lastCol, int maxRows);

template <typename T>
void PartsGenPermZeroDistinct(T* mat, const std::vector<T> &v,
                              std::vector<int> &z, int width,
                              int lastElem, int lastCol, int nRows);

void PartsDistinct(int* mat, std::vector<int> &z, int width,
                   int lastElem, int lastCol, int nRows);

void PartsDistinct(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                   int strt, int width, int lastElem,
                   int lastCol, int nRows);

void PartsPermDistinct(int* mat, std::vector<int> &z, int width,
                       int lastElem, int lastCol, int nRows);

void PartsPermZeroDistinct(int* mat, std::vector<int> &z, int width,
                           int lastElem, int lastCol, int nRows);

#endif