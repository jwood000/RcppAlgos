#ifndef GENERAL_PARTITIONS_H
#define GENERAL_PARTITIONS_H

#include "CombPermUtils.h"

int PartitionsMultiSet(int m, const std::vector<int64_t> &v, int64_t target, 
                       int lastElem, int lastCol, int maxRows, bool isComb,
                       const std::vector<int> &Reps, std::vector<int64_t> &partitionsVec);
    
int PartitionsRep(int m, const std::vector<int64_t> &v, int64_t target, int lastElem,
                  int lastCol, int maxRows, bool isComb, std::vector<int64_t> &partitionsVec);

int PartitionsDistinct(int m, const std::vector<int64_t> &v, int64_t target, int lastElem,
                       int lastCol, int maxRows, bool isComb, std::vector<int64_t> &partitionsVec);

#endif
