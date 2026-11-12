#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/Phi.h"
#include "CppConvert.h"
#include "cpp11/as.hpp"
#include <algorithm>
#include <cstdint>
#include <limits>

namespace {
    void preparePartialSieveInputs(
        SEXP Rn, SEXP Ra, SEXP RNumThreads, SEXP RmaxThreads,
        std::int64_t &n, std::int64_t &a, int &nThreads
    ) {

        nThreads = 1;
        double dblNum;
        double dblA;

        // Keep checkWhole enabled; decimalFraction permits zero in the sign check.
        CppConvert::convertPrimitive(Rn, dblNum, VecType::Numeric,
                                     "n", true, true, false, true);
        CppConvert::convertPrimitive(Ra, dblA, VecType::Numeric,
                                     "a", true, true, false, true);

        n = static_cast<std::int64_t>(dblNum);
        a = static_cast<std::int64_t>(dblA);

        int maxThreads = 1;
        CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                     VecType::Integer, "maxThreads");

        if (!Rf_isNull(RNumThreads)) {
            CppConvert::convertPrimitive(RNumThreads, nThreads,
                                         VecType::Integer, "nThreads");
        }

        nThreads = std::min(nThreads, maxThreads);
    }
}

[[cpp11::register]]
SEXP PrimePhiCpp(SEXP Rn, SEXP Ra, SEXP RNumThreads, SEXP RmaxThreads) {

    std::int64_t n;
    std::int64_t a;
    int nThreads;

    preparePartialSieveInputs(
        Rn, Ra, RNumThreads, RmaxThreads, n, a, nThreads
    );

    using namespace PrimeCounting;
    std::int64_t result;

    if (n == 0 || a == 0) {
        result = n;
    } else if (a >= n) {
        result = 1;
    } else if (a <= 6) {
        // The existing periodic evaluator needs only the first six primes.
        // In particular, Phi(n, 1) should not allocate a sqrt(n) lookup table.
        const PhiTables tables(13);
        result = Phi(n, a, tables);
    } else {
        const PhiTables tables(PartialSieveDetail::IntegerSqrt(n));
        const std::int64_t b = tables.primeCount();
        const std::int64_t survivors = PhiParallel(n, std::min(a, b),
                                                   nThreads, tables);

        // Once every prime through sqrt(n) has been excluded, only 1 and primes
        // remain. Each further excluded prime removes one survivor, until only
        // 1 remains. This covers a > b without storing primes through p[a].
        result = a <= b ? survivors :
            std::max<std::int64_t>(1, survivors - (a - b));
    }

    if (result > std::numeric_limits<int>::max()) {
        return cpp11::as_sexp(static_cast<double>(result));
    } else {
        return cpp11::as_sexp(static_cast<int>(result));
    }
}

[[cpp11::register]]
SEXP SemiPrimeCountCpp(SEXP Rn, SEXP Ra, SEXP RNumThreads, SEXP RmaxThreads) {

    std::int64_t n;
    std::int64_t a;
    int nThreads;

    preparePartialSieveInputs(
        Rn, Ra, RNumThreads, RmaxThreads, n, a, nThreads
    );

    using namespace PrimeCounting;
    const std::int64_t sqrtLimit = PartialSieveDetail::IntegerSqrt(n);

    // Since pi(sqrtLimit) <= sqrtLimit, this sufficient emptiness check
    // avoids preparing a table when every possible smaller factor is excluded.
    if (n < 4 || a >= sqrtLimit) {
        return cpp11::as_sexp(0);
    }

    const PhiTables tables(sqrtLimit);

    // The smaller factor of a nondecreasing prime pair is at most sqrt(n).
    if (a >= tables.primeCount()) {
        return cpp11::as_sexp(0);
    }

    const std::int64_t result = P2(n, a, nThreads, tables);

    if (result > std::numeric_limits<int>::max()) {
        return cpp11::as_sexp(static_cast<double>(result));
    } else {
        return cpp11::as_sexp(static_cast<int>(result));
    }
}

[[cpp11::register]]
SEXP TriPrimeCountCpp(SEXP Rn, SEXP Ra, SEXP RNumThreads, SEXP RmaxThreads) {

    std::int64_t n;
    std::int64_t a;
    int nThreads;

    preparePartialSieveInputs(
        Rn, Ra, RNumThreads, RmaxThreads, n, a, nThreads
    );

    using namespace PrimeCounting;
    const std::int64_t cubeRoot = PartialSieveDetail::IntegerCbrt(n);

    // Since pi(cubeRoot) <= cubeRoot, this sufficient emptiness check
    // avoids preparing a table when every possible smallest factor is excluded.
    if (n < 8 || a >= cubeRoot) {
        return cpp11::as_sexp(0);
    }

    std::int64_t result;
    std::int64_t tableLimit;
    {
        // The first factor of a nondecreasing triple is at most cbrt(n).
        // This small initial table determines emptiness and p[a + 1] before
        // attempting the potentially much larger table required by P3.
        const PhiTables tables(cubeRoot);

        if (a >= tables.primeCount()) {
            return cpp11::as_sexp(0);
        }

        const std::int64_t firstPrime = tables.primes()[a + 1];
        // Staged integer division gives floor(n / firstPrime^2) without
        // forming the square.
        tableLimit = std::max(cubeRoot, n / firstPrime / firstPrime);

        if (tableLimit == cubeRoot) {
            result = P3(n, a, nThreads, tables);
        }
    }

    if (tableLimit != cubeRoot) {
        // PhiTables supports coverage through floor(sqrt(2^53 - 1)). Unlike
        // Lehmer's cutoff, a small user-supplied a can require more than this.
        // Check the existing table limit before allocating the enlarged table.
        const std::int64_t maxTableLimit =
            PartialSieveDetail::IntegerSqrt(PartialSieveDetail::MAX_X);

        if (tableLimit > maxTableLimit) {
            cpp11::stop("triPrimeCount requires a lookup table beyond the supported "
                        "limit of %lld; reduce n or increase a",
                        static_cast<long long>(maxTableLimit));
        }

        const PhiTables tables(tableLimit);
        result = P3(n, a, nThreads, tables);
    }

    if (result > std::numeric_limits<int>::max()) {
        return cpp11::as_sexp(static_cast<double>(result));
    } else {
        return cpp11::as_sexp(static_cast<int>(result));
    }
}
