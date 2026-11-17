#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/Phi.h"
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "CppConvert/Constants.h"
#include <algorithm>
#include <stdexcept>

namespace PrimeCounting {

    std::int64_t PiMeissel(std::int64_t n, int nThreads, int maxThreads) {
        if (n < 0 || n > Significand53Int) {
            throw std::out_of_range(
                "PiMeissel n must be between 0 and 2^53 - 1"
            );
        }

        if (n < 2) {
            return 0;
        }

        const std::int64_t root = PartialSieveDetail::IntegerSqrt(n);
        const PhiTables tables(root);
        return PiMeissel(n, std::min(nThreads, maxThreads), tables);
    }

    std::int64_t PiMeissel(std::int64_t n, int nThreads,
                           const PhiTables& tables) {
        if (n < 0 || n > Significand53Int) {
            throw std::out_of_range(
                "PiMeissel n must be between 0 and 2^53 - 1"
            );
        }

        if (n < 2) {
            return 0;
        }

        const std::int64_t root = PartialSieveDetail::IntegerSqrt(n);

        if (tables.limit() < root) {
            throw std::out_of_range(
                "PiMeissel requires tables through floor(sqrt(n))"
            );
        }

        const std::int64_t a = tables.pi()[PartialSieveDetail::IntegerCbrt(n)];

        // Every surviving prime factor exceeds cbrt(n), so a composite
        // survivor has exactly two factors, counted with multiplicity:
        // Phi(n, a) = 1 + (pi(n) - a) + P2(n, a).
        // Complete each phase before starting the next, sharing the prepared
        // tables and worker budget without nesting parallel launches.
        const std::int64_t phi = PhiParallel(n, a, nThreads, tables);
        const std::int64_t p2 = P2(n, a, nThreads, tables);
        return phi + a - 1 - p2;
    }
}
