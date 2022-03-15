#ifndef PARTITIONS_REP_H
#define PARTITIONS_REP_H

#include "RMatrix.h"
#include <vector>

template <typename T>
void PartsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                 int width, int lastElem, int lastCol, int nRows);

template <typename T>
void PartsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                 std::vector<int> &z, int strt, int width, int lastElem,
                 int lastCol, int nRows);

template <typename T>
void PartsGenPermRep(T*, const std::vector<T> &v, std::vector<int> &z,
                     int width, int lastElem, int lastCol, int nRows);

template <typename T>
void PartsGenRep(std::vector<T> &partsVec, const std::vector<T> &v,
                 std::vector<int> &z, int width, int nRows, bool IsComb);

void PartsRep(int* mat, std::vector<int> &z, int width,
              int lastElem, int lastCol, int nRows);

void PartsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
              int strt, int width, int lastElem, int lastCol, int nRows);

void PartsPermRep(int* mat, std::vector<int> &z, int width,
                  int lastElem, int lastCol, int nRows);

#endif
