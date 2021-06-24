#ifndef PARTITIONS_DISTINCT_H
#define PARTITIONS_DISTINCT_H

#include <vector>

template <typename T>
void PartsGenDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, int width, int lastElem,
                      int lastCol, int strt, int nRows);

template <typename T>
void PartsGenPermDistinct(T* mat, const std::vector<T> &v,
                          std::vector<int> &z, int width,
                          int lastElem, int lastCol, int maxRows);

template <typename T>
void PartsGenPermZeroDistinct(T* mat, const std::vector<T> &v,
                              std::vector<int> &z, int width,
                              int lastElem, int lastCol, int nRows);

void PartsDistinct(int* mat, std::vector<int> &z, int width, int boundary,
                   int lastCol, int edge, int strt, int nRows);

void PartsPermDistinct(int* mat, std::vector<int> &z, int width,
                       int boundary, int lastCol, int edge, int nRows);

void PartsPermZeroDistinct(int* mat, std::vector<int> &z, int width,
                           int boundary, int lastCol, int edge, int nRows);

#endif