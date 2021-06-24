#ifndef PARTITIONS_REP_H
#define PARTITIONS_REP_H

#include <vector>

template <typename T>
void PartsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                 int width, int lastElem, int lastCol, int strt, int nRows);

template <typename T>
void PartsGenPermRep(T*, const std::vector<T> &v, std::vector<int> &z,
                     int width, int lastElem, int lastCol, int nRows);

template <typename T>
void PartsGenPermRep(std::vector<T> &partsVec, const std::vector<T> &v,
                     std::vector<int> &z, int width, int nRows);

void PartsRep(int* mat, std::vector<int> &z, int width, int boundary,
              int edge, int lastCol, int strt, int nRows);

void PartsPermRep(int* mat, std::vector<int> &z, int width,
                  int boundary, int edge, int lastCol, int nRows);

#endif
