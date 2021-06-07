#ifndef PARTITIONS_DISTINCT_H
#define PARTITIONS_DISTINCT_H

#include <vector>

template <typename T>
void PartsGenDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, int width, int lastElem,
                      int lastCol, int strt, int nRows);

template <typename T>
void PartsGenPermDistinct(std::vector<T> &partitionsVec,
                          const std::vector<T> &v, std::vector<int> &z,
                          int width, int lastElem, int lastCol, int maxRows);

void PartsDistinct(int* mat, std::vector<int> &z, int width, int boundary,
                   int lastCol, int edge, int strt, int nRows);

// mIsNull && IncludeZero
void PartsPermDistinct(int* mat, std::vector<int> &z, int width,
                       int boundary, int lastCol, int edge, int nRows);

// !mIsNull || !IncludeZero
void PartsLenPermDistinct(int* mat, std::vector<int> &z,
                          int width, int boundary, int lastCol,
                          int edge, int nRows);

#endif