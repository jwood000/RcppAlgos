#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/Phi.h"
#include <algorithm>
#include <cmath>

// "PiLegendre" based off of the highly optimized
// "pi_legendre.cpp" algorithm by Kim Walisch, which calculates
// the numbers of primes less than n using Legendre's formula.
// Kim Walisch's official github repo for "pi_legendre" is:
//                      https://github.com/kimwalisch/primecount

namespace PrimeCounting {

    // Historical timings before shared-Phi integration (nThreads = 8).
    // Counts verified by Kim Walisch's primecount library.
    //  10^9 -->>          50,847,534   -->>   192.8 microseconds
    // 10^10 -->>         455,052,511   -->>   800.9 microseconds
    // 10^11 -->>       4,118,054,813   -->>   2.799 milliseconds
    // 10^12 -->>      37,607,912,018   -->>   14.72 milliseconds
    // 10^13 -->>     346,065,536,839   -->>   104.3 milliseconds
    // 10^14 -->>   3,204,941,750,802   -->>   0.770 seconds
    // 10^15 -->>  29,844,570,422,669   -->>   5.779 seconds
    // MAX VALUE (2^53 - 1) -->>
    //            252,252,704,148,404   -->>  42.165 seconds

    std::int64_t PiLegendre(std::int64_t n, int nThreads, int maxThreads) {

        const std::int64_t sqrtBound = std::sqrt(static_cast<double>(n));
        const PhiTables tables(sqrtBound);

        nThreads = std::min(nThreads, maxThreads);

        const std::int64_t piSqrt = tables.primeCount();
        const std::int64_t phiSqrt = PhiParallel(n, piSqrt, nThreads, tables);
        const std::int64_t int64result = piSqrt + phiSqrt - 1;

        return int64result;
    }
}
