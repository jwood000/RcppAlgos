#ifndef NTH_PARTITION_H
#define NTH_PARTITION_H

#include "Partitions/PartitionsTypes.h"
#include <gmp.h>

using nthPartsPtr = std::vector<int> (*const)(int n, int m, int cap, int k,
                                              double dblIdx, mpz_t mpzIdx);

nthPartsPtr GetNthPartsFunc(PartitionType ptype, bool IsGmp, bool IsComp);

#endif
