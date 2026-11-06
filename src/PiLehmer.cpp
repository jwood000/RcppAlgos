#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/Phi.h"
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include <algorithm>
#include <stdexcept>

namespace PrimeCounting {

    std::int64_t PiLehmer(std::int64_t n, int nThreads, int maxThreads) {
        if (n < 0 || n > PartialSieveDetail::MAX_X) {
            throw std::out_of_range("PiLehmer n must be between 0 and 2^53 - 1");
        }
        if (n < 2) return 0;

        const std::int64_t root = PartialSieveDetail::IntegerSqrt(n);
        const PhiTables tables(root);
        // floor(sqrt(floor(sqrt(n)))) is exactly floor(n^(1/4)).
        const std::int64_t a = tables.pi()[PartialSieveDetail::IntegerSqrt(root)];
        nThreads = std::min(nThreads, maxThreads);

        // Surviving prime factors exceed n^(1/4), so there are at most three:
        // Phi(n, a) = 1 + (pi(n) - a) + P2(n, a) + P3(n, a).
        // The same cutoff puts every P3 pi lookup below sqrt(n), within the
        // shared table. Sequential phases keep their worker groups separate.
        const std::int64_t phi = PhiParallel(n, a, nThreads, tables);
        const std::int64_t p2 = P2(n, a, nThreads, tables);
        const std::int64_t p3 = P3(n, a, nThreads, tables);
        return phi + a - 1 - p2 - p3;
    }
}
