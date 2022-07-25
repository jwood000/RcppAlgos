#ifndef RANK_PARTITION_H
#define RANK_PARTITION_H

#include "Partitions/PartitionsTypes.h"
#include <gmp.h>

using rankPartsPtr = void (*const)(std::vector<int>::iterator iter,
                           int n, int m, int cap, int k,
                           double &dblIdx, mpz_t mpzIdx);

rankPartsPtr GetRankPartsFunc(PartitionType ptype, bool IsGmp, bool IsComp);

#endif
