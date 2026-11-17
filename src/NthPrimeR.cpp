#include "NumbersUtils/NthPrime.h"
#include "CppConvert.h"
#include "cpp11/as.hpp"
#include <algorithm>
#include <limits>

[[cpp11::register]]
SEXP NthPrimeCpp(SEXP Rn, SEXP RNumThreads, SEXP RmaxThreads, SEXP RmaxCores) {

    double dblNum;
    CppConvert::convertPrimitive(Rn, dblNum, VecType::Numeric, "n");

    const std::int64_t n = static_cast<std::int64_t>(dblNum);

    if (n > PrimeCounting::MAX_NTH_PRIME_INDEX) {
        cpp11::stop(
            "n must be at most 252252704148404 so the nth prime is at most 2^53 - 1"
        );
    }

    int nThreads = 1;
    int maxThreads = 1;
    int maxCores = 1;

    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");
    CppConvert::convertPrimitive(RmaxCores, maxCores,
                                 VecType::Integer, "maxCores");

    if (!Rf_isNull(RNumThreads)) {
        CppConvert::convertPrimitive(RNumThreads, nThreads,
                                     VecType::Integer, "nThreads");
    }

    nThreads = std::min(nThreads, maxThreads);
    const std::int64_t result = PrimeCounting::NthPrime(n, nThreads,
                                                      maxThreads, maxCores);

    if (result > std::numeric_limits<int>::max()) {
        return cpp11::as_sexp(static_cast<double>(result));
    } else {
        return cpp11::as_sexp(static_cast<int>(result));
    }
}
