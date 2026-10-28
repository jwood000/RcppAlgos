#pragma once

#include "Partitions/PartitionsTypes.h"
#include <gmpxx.h>
#include <vector>

using rankPartsFunc = void(
    std::vector<int>::iterator iter, int n, int m, int cap, int k,
    double &dblIdx, mpz_class &mpzIdx
);
using rankPartsPtr = rankPartsFunc *const;

rankPartsPtr GetRankPartsFunc(PartitionType ptype, bool IsGmp);

// Partition implementations, with double/GMP counterparts paired.
rankPartsFunc rankPartsRepLen;
rankPartsFunc rankPartsRepLenGmp;
rankPartsFunc rankPartsRepShort;
rankPartsFunc rankPartsRepShortGmp;
rankPartsFunc rankPartsRep;
rankPartsFunc rankPartsRepGmp;
rankPartsFunc rankPartsRepCap;
rankPartsFunc rankPartsRepCapGmp;
rankPartsFunc rankPartsDistinctLen;
rankPartsFunc rankPartsDistinctLenGmp;
rankPartsFunc rankPartsDistinctOneZero;
rankPartsFunc rankPartsDistinctOneZeroGmp;
rankPartsFunc rankPartsDistinctMZ;
rankPartsFunc rankPartsDistinctMZGmp;
rankPartsFunc rankPartsDistinctCap;
rankPartsFunc rankPartsDistinctCapGmp;
rankPartsFunc rankPartsDistinctCapMZ;
rankPartsFunc rankPartsDistinctCapMZGmp;
