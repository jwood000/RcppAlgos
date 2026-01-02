#pragma once

#include "Partitions/PartitionsTypes.h"
#include <gmpxx.h>

// The variable k is strtLen
using nthPartsPtr = std::vector<int> (*const)(int n, int m, int cap, int k,
                                      double dblIdx, const mpz_class &mpzIdx);

nthPartsPtr GetNthPartsFunc(PartitionType ptype, bool IsGmp);
nthPartsPtr GetNthPartsFuncOrStop(PartitionType ptype, bool IsGmp);
