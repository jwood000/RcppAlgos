#pragma once

#include "Partitions/PartitionsTypes.h"
#include "RMatrix.h"

int PartsStdManager(
    int* mat, std::vector<int> &z, int width, int lastElem,
    int lastCol, int nRows, PartitionType ptype
);

int PartsStdParallel(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                     int strt, int width, int lastElem, int lastCol,
                     int nRows, PartitionType ptype);

template <typename T>
int PartsGenParallel(RcppParallel::RMatrix<T> &mat,
                     const std::vector<T> &v, std::vector<int> &z, int strt,
                     int width, int lastElem, int lastCol, int nRows,
                     PartitionType ptype);

template <typename T>
int PartsGenManager(T* mat, const std::vector<T> &v, std::vector<int> &z,
                    int width, int lastElem, int lastCol, int nRows,
                    PartitionType ptype);

template <typename T>
int PartsGenManager(std::vector<T> &partsVec, const std::vector<T> &v,
                    const std::vector<int> &Reps, std::vector<int> &z,
                    int width, int nRows, PartitionType ptype);
