#ifndef PARTITIONS_GENERAL_H
#define PARTITIONS_GENERAL_H

#include "CleanConvert.h"
#include <vector>
#include <cstdint>

template <typename typeVector>
void PartitionsMultiSet(int m, const std::vector<typeVector> &v, std::vector<int> &z,
                        int lastElem, int lastCol, int maxRows, bool isComb,
                        const std::vector<int> &Reps, std::vector<typeVector> &partitionsVec);

template <typename typeRcpp, typename typeVector>
void PartitionsRep(int m, const std::vector<typeVector> &v, std::vector<int> &z,
                   int lastElem, int lastCol, int maxRows, bool isComb,
                   std::vector<typeVector> &partitionsVec, typeRcpp &matRcpp);

template <typename typeRcpp, typename typeVector>
void PartitionsDistinct(int m, const std::vector<typeVector> &v, std::vector<int> &z,
                        int lastElem, int lastCol, int maxRows, bool isComb,
                        std::vector<typeVector> &partitionsVec, typeRcpp &matRcpp);

#endif
