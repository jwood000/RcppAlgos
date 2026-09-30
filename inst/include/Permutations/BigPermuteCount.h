#pragma once

#include <vector>
#include <gmpxx.h>

// `v` must be sorted in non-decreasing order. When includeZero is false,
// only leading zeros are excluded from the count.
void NumPermsWithRepGmp(
    mpz_class &result, const std::vector<int> &v, bool includeZero = true
);
void NumPermsNoRepGmp(mpz_class &result, int n, int m);
void MultisetPermRowNumGmp(mpz_class &result, int n, int m,
                           const std::vector<int> &myReps);
