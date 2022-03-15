#ifndef THREAD_SAFE_PARTS_H
#define THREAD_SAFE_PARTS_H

#include "Partitions/PartitionsTypes.h"
#include <gmp.h>

void StandardPartitions(int* mat, std::vector<int> &z, PartitionType ptype,
                        double lower, mpz_t lowerMpz, int nCols,
                        int width, int nRows, int nThreads, int lastCol,
                        int lastElem, int tar, int strtLen, int cap,
                        bool IsRep, bool IsMult, bool IsGmp,
                        bool IsComb, bool includeZero);

template <typename T>
void GeneralPartitions(T* mat, const std::vector<T> &v, std::vector<int> &z,
                       const PartDesign &part, double lower, mpz_t lowerMpz,
                       int nCols, int nRows, int nThreads, int lastCol,
                       int lastElem, int strtLen, int cap, bool IsComb);

#endif