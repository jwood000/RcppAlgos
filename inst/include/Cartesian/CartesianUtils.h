#pragma once

#include "cpp11/strings.hpp"

#include "SetUpUtils.h"
#include <numeric>

double CartesianCount(const std::vector<int> &lenGrps);

void CartesianCountGmp(mpz_class &result, const std::vector<int> &lenGrps);

std::vector<int> nthProduct(double dblIdx, const std::vector<int> &lenGrp);

std::vector<int> nthProductGmp(const mpz_class &mpzIdx,
                               const std::vector<int> &lenGrp);

bool nextProduct(const std::vector<int> &lenGrps,
                 std::vector<int> &z, int m);

bool prevProduct(const std::vector<int> &lenGrps,
                 std::vector<int> &z, int m);

void GetStartProd(
    const std::vector<int> &lenNxtPr, std::vector<int> &z,
    mpz_class &lowerMpz, double &lower, int stepSize, bool IsGmp
);

void CartesianInitialPrep(
    cpp11::list RList, std::vector<int> &IsFactor,
    std::vector<int> &lenGrps, int nCols
);

void ProductPrepare(
    cpp11::list RList, const std::vector<int> &IsFactor,
    const std::vector<int> &lenGrps, std::vector<std::vector<int>> &myVec,
    cpp11::writable::strings &charVec, std::vector<Rcomplex> &cmplxVec,
    std::vector<Rbyte> &rawVec, std::vector<double> &dblVec,
    std::vector<int> &intVec, std::vector<int> &boolVec,
    std::vector<int> &typeCheck, VecType &myType, int nCols, bool &IsDF
);
