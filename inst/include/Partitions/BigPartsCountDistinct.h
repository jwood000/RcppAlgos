#pragma once

#include <gmpxx.h>
#include <vector>

void CountPartsDistLenRstrctd(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen = 0
);

void CountPartsDistinctLen(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);

void CountPartsDistinct(
    mpz_class &res, int n, int m,
    const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);

void CountPartsDistinctMultiZero(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed, int strtLen
);

void CountPartsDistinctRstrctdMZ(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen
);

void CountCompDistLenRstrctd(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen = 0
);

void CountPartsPermDistinctRstrctdMZ(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen
);

void CountCompsDistinctLen(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);

void CountCompsDistinctMultiZero(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed, int strtLen
);

void CountCompsDistinctMZWeak(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed, int strtLen
);

void CountCompsDistinctRstrctdMZ(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen
);
