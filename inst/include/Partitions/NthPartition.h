#pragma once

#include "Partitions/PartitionsTypes.h"
#include <gmpxx.h>
#include <vector>

// The variable k is strtLen. Reps contains the input multiplicities.
using nthPartsFunc = std::vector<int>(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> &Reps
);
using nthPartsPtr = nthPartsFunc *const;

// Unranking uses zero-based indices and skips whole lexicographic blocks.
// A zero completion count denotes an empty block and must be skipped normally.
// Callers supply a positive width and validate the index against the total
// count before unranking.
//
// Positive/mapped algorithms return offset values: res[i] = part_i - 1.
// Their final element is (target - width) - sum(res). Zero-based algorithms
// instead use target - sum(res). Both families share these conventions.

nthPartsPtr GetNthPartsFunc(PartitionType ptype, bool IsGmp);
nthPartsPtr GetNthPartsFuncOrStop(PartitionType ptype, bool IsGmp);

// Partition implementations, with double/GMP counterparts paired.
nthPartsFunc nthPartsRepLen;
nthPartsFunc nthPartsRepLenGmp;
nthPartsFunc nthPartsRepShort;
nthPartsFunc nthPartsRepShortGmp;
nthPartsFunc nthPartsRep;
nthPartsFunc nthPartsRepGmp;
nthPartsFunc nthPartsRepCap;
nthPartsFunc nthPartsRepCapGmp;
nthPartsFunc nthPartsDistinctLen;
nthPartsFunc nthPartsDistinctLenGmp;
nthPartsFunc nthPartsDistinctMZ;
nthPartsFunc nthPartsDistinctMZGmp;
nthPartsFunc nthPartsDistinctCap;
nthPartsFunc nthPartsDistinctCapGmp;
nthPartsFunc nthPartsDistinctCapMZ;
nthPartsFunc nthPartsDistinctCapMZGmp;

// Multiset contract (double and GMP): m > 0, cap == Reps.size() > 0,
// and every input multiplicity is positive. n is the target for normalized
// values 1..cap; the returned values are zero-based bucket indices.
// The active index must satisfy 0 <= index < the count for (n, m, Reps).
// Public sampleVec indices are one-based and validated before conversion;
// direct C++ callers must establish these preconditions themselves.
nthPartsFunc nthPartsMultiset;
nthPartsFunc nthPartsMultisetGmp;
