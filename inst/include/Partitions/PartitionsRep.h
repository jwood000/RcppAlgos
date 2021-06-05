#ifndef PARTITIONS_REP_H
#define PARTITIONS_REP_H

#include <vector>

template <typename T>
void PartsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                 int m, int lastElem, int lastCol, int strt, int nRows);

template <typename T>
void PartsGenPermRep(std::vector<T> &partitionsVec,
                     const std::vector<T> &v, std::vector<int> &z,
                     int m, int lastElem, int lastCol, int maxRows);

void PartsRep(int* mat, std::vector<int> &z, int m, int boundary,
              int edge, int lastCol, int strt, int nRows);

void PartsPermRep(int* mat, std::vector<int> &z, int m,
                  int boundary, int edge, int lastCol, int nRows);

#endif
