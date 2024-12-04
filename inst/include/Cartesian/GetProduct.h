#pragma once

#include "cpp11/logicals.hpp"
#include "cpp11/matrix.hpp"

#include "RMatrix.h"
#include <functional>
#include <algorithm>
#include <cstdint>
#include <thread>

#include "Cartesian/CartesianUtils.h"

SEXP GetProduct(
    const std::vector<int> &idx, const std::vector<int> &typeCheck,
    const std::vector<int> &IsFactor, const cpp11::list &RList,
    const std::vector<int> &intVec, const std::vector<double> &dblVec,
    const std::vector<int> &boolVec, const std::vector<Rcomplex> &cmplxVec,
    const std::vector<Rbyte> &rawVec, const cpp11::strings &charVec,
    const std::vector<int> &lenGrps, std::vector<int> &z,
    const std::vector<double> &mySamp, const std::vector<mpz_class> &myBigSamp,
    double lower, mpz_class &lowerMpz, int nRows, int nCols, bool IsDF,
    int nThreads, bool Parallel, bool IsGmp, bool IsSample
);
