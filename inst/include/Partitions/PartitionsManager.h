#ifndef PARTITIONS_MANAGER_H
#define PARTITIONS_MANAGER_H

#include "Partitions/PartitionsTypes.h"

void PartsStdManager(int* mat, std::vector<int> &z, int width,
                     int lastCol, int boundary, int edge,
                     int strt, int nRows, bool IsComb, bool IsRep);

template <typename T>
void PartsGenManager(T* mat, const PartDesign &part, const std::vector<T> &v,
                     std::vector<int> &z, int width, int lastCol,
                     int lastElem, int boundary, int edge, int strt,
                     int nRows, bool IsComb, bool IsRep);

template <typename T>
void PartsGenManager(std::vector<T> &partsVec, const std::vector<T> &v,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     PartitionType ptype, int width, int nRows, bool IsComb);

#endif