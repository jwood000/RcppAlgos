#include "NumbersUtils/PrimeCount.h"
#include "CppConvert.h"
#include "cpp11/as.hpp"
#include <limits>
#include <string>

[[cpp11::register]]
SEXP PrimeCountCpp(SEXP Rn, SEXP RNumThreads, SEXP RmaxThreads,
                   std::string method) {
    double dblNum;
    CppConvert::convertPrimitive(Rn, dblNum, VecType::Numeric, "n");
    const std::int64_t n = static_cast<std::int64_t>(dblNum);

    int nThreads = 1;
    int maxThreads = 1;
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

    using Counter = std::int64_t (*)(std::int64_t, int, int);
    Counter count = nullptr;

    // Auto initially retains the historical sieve/Legendre policy. Explicit
    // choices use the named algorithm for every accepted n, including small n.
    if (method == "auto" || method == "legendre") {
        count = PrimeCounting::PiLegendre;
    } else if (method == "meissel") {
        count = PrimeCounting::PiMeissel;
    } else if (method == "lehmer") {
        count = PrimeCounting::PiLehmer;
    } else {
        cpp11::stop("Unknown prime-counting method");
    }

    if (method == "auto" && n < 100000) {
        if (n < 10) {
            if (n == 1)
                return cpp11::as_sexp(0);
            else if (n == 2)
                return cpp11::as_sexp(1);
            else if (n < 5)
                return cpp11::as_sexp(2);
            else if (n < 7)
                return cpp11::as_sexp(3);
            else
                return cpp11::as_sexp(4);
        }

        return cpp11::as_sexp(static_cast<int>(PrimeCounting::PiPrime(n)));
    }

    // Preserve the existing small-input contract: ignore nThreads below
    // 100000. Explicit algorithms also run serially for these inputs.
    if (n >= 100000 && !Rf_isNull(RNumThreads)) {
        CppConvert::convertPrimitive(RNumThreads, nThreads,
                                       VecType::Integer, "nThreads");
    }

    const std::int64_t result = count(n, nThreads, maxThreads);

    if (result > std::numeric_limits<int>::max()) {
        return cpp11::as_sexp(static_cast<double>(result));
    } else {
        return cpp11::as_sexp(static_cast<int>(result));
    }
}
