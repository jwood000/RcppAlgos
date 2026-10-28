#pragma once

#include "Partitions/PartitionsTypes.h"
#include <gmpxx.h>
#include <vector>

using rankPartsFunc = void(
    std::vector<int>::iterator iter, int n, int m, int cap, int k,
    double &dblIdx, mpz_class &mpzIdx, const std::vector<int> &Reps
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

// Multiset contract (double and GMP): m > 0, cap == Reps.size() > 0,
// and every input multiplicity is positive. iter addresses m nondecreasing
// bucket indices in [0, cap), summing to n - m, with no bucket used more
// often than Reps permits. The result is zero-based.
// GetRankPart checks membership/target; SetUpRank sorts and checks usage.
// Direct C++ callers must establish these preconditions themselves.
rankPartsFunc rankPartsMultiset;
rankPartsFunc rankPartsMultisetGmp;
