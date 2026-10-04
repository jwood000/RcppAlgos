#pragma once

#include <gmpxx.h>
#include <vector>

void CountPartsMultiset(
    mpz_class &res, int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen = 0
);

void CountCompsMultiset(
    mpz_class &res, int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen = 0
);

void CountCompsMultisetZero(
    mpz_class &res, int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
);

void CountCompsMultisetWeak(
    mpz_class &res, int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
);
