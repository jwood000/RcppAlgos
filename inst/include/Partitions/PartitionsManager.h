#ifndef PARTITIONS_MANAGER_H
#define PARTITIONS_MANAGER_H

#include "Partitions/PartitionsTypes.h"

template <typename T>
void PartsGenManager(T* mat, const std::vector<T> &v,
                     std::vector<T> &partVec,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     int width, int lastCol, int lastElem, int boundary,
                     int edge, int strt, int nRows, bool IsComb,
                     bool IsRep, bool IsMult);

void PartsStdManager(int* mat, std::vector<int> &z, int width,
                     int lastCol, int boundary, int edge,
                     int strt, int nRows, bool IsComb, bool IsRep);

#endif