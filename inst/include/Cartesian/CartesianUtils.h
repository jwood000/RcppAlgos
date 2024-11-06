#pragma once

#include <cstdlib>
#include <numeric>
#include <vector>
#include <gmpxx.h>

double CartesianCount(const std::vector<int> &lenGrps);

void CartesianCountGmp(mpz_class &result, const std::vector<int> &lenGrps);

std::vector<int> nthProduct(double dblIdx, const std::vector<int> &lenGrp);

std::vector<int> nthProductGmp(const mpz_class &mpzIdx,
                               const std::vector<int> &lenGrp);

bool nextProduct(const std::vector<int> &lenGrps,
                 std::vector<int> &z, int m);

void GetStartProd(
    const std::vector<int> &lenNxtPr, std::vector<int> &z,
    mpz_class &lowerMpz, double &lower, int stepSize, bool IsGmp
);
