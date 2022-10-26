#pragma once

#include "Partitions/PartitionsTypes.h"
#include "RMatrix.h"

void PartsStdManager(int* mat, std::vector<int> &z, int width,
                     int lastElem, int lastCol, int nRows, bool IsComb,
                     bool IsRep, bool IsComp, bool zero_spesh);

void PartsStdParallel(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                      int strt, int width, int lastElem, int lastCol,
                      int nRows, bool IsRep, bool IsComp, bool zero_spesh);

template <typename T>
void PartsGenParallel(RcppParallel::RMatrix<T> &mat,
                      const std::vector<T> &v, std::vector<int> &z,
                      int strt, int width, int lastElem, int lastCol,
                      int nRows, bool IsRep, bool IsComp, bool zero_spesh);

template <typename T>
void PartsGenManager(T* mat, const std::vector<T> &v, std::vector<int> &z,
                     int width, int lastElem, int lastCol, int nRows,
                     bool IsComb, bool IsRep, bool IsComp, bool zero_spesh);

template <typename T>
void PartsGenManager(std::vector<T> &partsVec, const std::vector<T> &v,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     PartitionType ptype, int width, int nRows, bool IsComb);
