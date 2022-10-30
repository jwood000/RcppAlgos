#pragma once

#include "Partitions/PartitionsTypes.h"

using rankPartsPtr = void (*const)(std::vector<int>::iterator iter,
                           int n, int m, int cap, int k,
                           double &dblIdx, mpz_class &mpzIdx);

rankPartsPtr GetRankPartsFunc(PartitionType ptype, bool IsGmp, bool IsComp);
