#ifndef NTH_PARTITION_H
#define NTH_PARTITION_H

#include <vector>
#include <gmp.h>

using nthPartsPtr = std::vector<int> (*const)(int n, int r,
                                      double dblIdx, mpz_t mpzIdx);

nthPartsPtr GetNthPartsFunc(bool IsMult, bool IsRep, bool IsGmp);

#endif
