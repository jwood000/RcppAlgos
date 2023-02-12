#pragma once

#include "Partitions/PartitionsTypes.h"
#include "CppConvert/GmpxxCopy.h"

using nthPartsPtr = std::vector<int> (*const)(int n, int m, int cap, int k,
                                      double dblIdx, const mpz_class &mpzIdx);

nthPartsPtr GetNthPartsFunc(PartitionType ptype, bool IsGmp, bool IsComp);
