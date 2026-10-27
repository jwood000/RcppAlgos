#pragma once

#include <gmpxx.h>
#include <vector>

// Each buffer must have at least (m + 1) * (n + 1) elements.
// Workers reset and reuse the buffers, including swapping their contents.
void CountPartsMultiset(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen = 0
);

void CountCompsMultiset(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen = 0
);

void CountCompsMultisetZero(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
);

void CountCompsMultisetWeak(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
);
