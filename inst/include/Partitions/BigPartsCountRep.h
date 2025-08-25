#pragma once

#include <gmpxx.h>
#include <vector>

void CountPartsRepLenRstrctd(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen = 0
);

void CountPartsRepLen(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);

void CountPartsRep(
    mpz_class &res, int n, int m,
    const std::vector<int> &allowed = std::vector<int>(), int strtLen = 0
);

void CountCompsRepLen(
    mpz_class &res, int n, int m,
    const std::vector<int> &allowed = std::vector<int>(), int strtLen = 0
);

void CountCompsRepZNotWk(
    mpz_class &res, int n, int m,
    const std::vector<int> &allowed = std::vector<int>(), int strtLen = 0
);
