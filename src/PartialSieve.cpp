#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/Phi.h"
#include "CppConvert.h"
#include "cpp11/as.hpp"
#include <algorithm>
#include <cstdint>
#include <limits>

namespace {

    void checkScalarNumeric(SEXP input, const char* name) {
        if (TYPEOF(input) != REALSXP && TYPEOF(input) != INTSXP) {
            cpp11::stop("%s must be of type numeric or integer", name);
        }
        if (Rf_xlength(input) != 1) {
            cpp11::stop("%s must be of length 1", name);
        }
    }

    std::int64_t getNonnegativeWhole(SEXP input, const char* name) {
        checkScalarNumeric(input, name);
        double value;
        // The shared converter checks finiteness, integrality and the 53-bit
        // range before conversion. Permit zero, then reject negative values.
        CppConvert::convertPrimitive(input, value, VecType::Numeric,
                                     name, true, true, true);
        if (value < 0) {
            cpp11::stop("%s must be a nonnegative whole number", name);
        }
        return static_cast<std::int64_t>(value);
    }

    int getThreadCount(SEXP requested, SEXP maximum) {
        int nThreads = 1;
        int maxThreads = 1;
        CppConvert::convertPrimitive(maximum, maxThreads,
                                     VecType::Integer, "maxThreads");
        if (!Rf_isNull(requested)) {
            checkScalarNumeric(requested, "nThreads");
            CppConvert::convertPrimitive(requested, nThreads,
                                         VecType::Integer, "nThreads");
        }
        return std::min(nThreads, maxThreads);
    }

    SEXP countResult(std::int64_t count) {
        if (count > std::numeric_limits<int>::max()) {
            return cpp11::as_sexp(static_cast<double>(count));
        }
        return cpp11::as_sexp(static_cast<int>(count));
    }

} // namespace

[[cpp11::register]]
SEXP PrimePhiCpp(SEXP Rn, SEXP Ra, SEXP RNumThreads, SEXP RmaxThreads) {
    const std::int64_t n = getNonnegativeWhole(Rn, "n");
    const std::int64_t a = getNonnegativeWhole(Ra, "a");
    const int nThreads = getThreadCount(RNumThreads, RmaxThreads);

    if (n == 0 || a == 0) return countResult(n);
    if (a >= n) return countResult(1);

    using namespace PrimeCounting;
    if (a <= 6) {
        // The existing periodic evaluator needs only the first six primes.
        // In particular, Phi(n, 1) should not allocate a sqrt(n) lookup table.
        const PhiTables tables(13);
        return countResult(Phi(n, a, tables));
    }

    const PhiTables tables(PartialSieveDetail::IntegerSqrt(n));
    const std::int64_t b = tables.primeCount();
    const std::int64_t survivors = PhiParallel(n, std::min(a, b),
                                               nThreads, tables);

    // Once every prime through sqrt(n) has been excluded, only 1 and primes
    // remain. Each further excluded prime removes one survivor, until only
    // 1 remains. This covers a > b without storing primes through p[a].
    return countResult(a <= b ? survivors :
        std::max<std::int64_t>(1, survivors - (a - b)));
}

[[cpp11::register]]
SEXP SemiPrimeCountCpp(SEXP Rn, SEXP Ra, SEXP RNumThreads, SEXP RmaxThreads) {
    const std::int64_t n = getNonnegativeWhole(Rn, "n");
    const std::int64_t a = getNonnegativeWhole(Ra, "a");
    const int nThreads = getThreadCount(RNumThreads, RmaxThreads);

    using namespace PrimeCounting;
    const std::int64_t sqrtLimit = PartialSieveDetail::IntegerSqrt(n);
    // Since pi(sqrtLimit) <= sqrtLimit, this sufficient emptiness check
    // avoids preparing a table when every possible smaller factor is excluded.
    if (n < 4 || a >= sqrtLimit) return countResult(0);

    const PhiTables tables(sqrtLimit);
    // The smaller factor of a nondecreasing prime pair is at most sqrt(n).
    if (a >= tables.primeCount()) return countResult(0);
    return countResult(P2(n, a, nThreads, tables));
}

[[cpp11::register]]
SEXP TriPrimeCountCpp(SEXP Rn, SEXP Ra, SEXP RNumThreads, SEXP RmaxThreads) {
    const std::int64_t n = getNonnegativeWhole(Rn, "n");
    const std::int64_t a = getNonnegativeWhole(Ra, "a");
    const int nThreads = getThreadCount(RNumThreads, RmaxThreads);

    using namespace PrimeCounting;
    const std::int64_t cubeRoot = PartialSieveDetail::IntegerCbrt(n);
    // Since pi(cubeRoot) <= cubeRoot, this sufficient emptiness check
    // avoids preparing a table when every possible smallest factor is excluded.
    if (n < 8 || a >= cubeRoot) return countResult(0);

    std::int64_t tableLimit;
    {
        // The first factor of a nondecreasing triple is at most cbrt(n).
        // This small initial table determines emptiness and p[a + 1] before
        // attempting the potentially much larger table required by P3.
        const PhiTables tables(cubeRoot);
        if (a >= tables.primeCount()) return countResult(0);
        const std::int64_t firstPrime = tables.primes()[a + 1];
        // Staged integer division gives floor(n / firstPrime^2) without
        // forming the square.
        tableLimit = std::max(cubeRoot, n / firstPrime / firstPrime);
        if (tableLimit == cubeRoot) {
            return countResult(P3(n, a, nThreads, tables));
        }
    }

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
    return countResult(P3(n, a, nThreads, tables));
}
