#pragma once

#include <cstdint>

namespace PrimeCounting {
    class PhiTables;

    // Count n <= x with exactly two/three prime factors, with multiplicity,
    // all greater than p[a]. a == 0 permits every prime. Both require
    // 0 <= x <= 2^53 - 1 and 0 <= a <= tables.primeCount(); invalid inputs or
    // insufficient table coverage throw std::out_of_range. Tables are borrowed
    // immutably. nThreads <= 1 requests serial execution; small jobs are serial.
    // Use at the chosen outer parallel level; worker computations stay serial.

    // For x >= 4, prepare tables through floor(sqrt(x)). Larger pi queries
    // use serial Phi with those same tables and its bounded shared cache.
    // No dense table through x / p[a + 1] or additional sieve is allocated.
    std::int64_t P2(std::int64_t x, std::int64_t a, int nThreads,
                    const PhiTables& tables);

    // For x >= 8, prepare tables through floor(cbrt(x)). If any triples exist,
    // also require coverage through floor(x / p[a + 1] / p[a + 1]). That bound
    // covers every inner prime and pi lookup. A sqrt(x) table suffices at
    // Lehmer's fourth-root cutoff, but not necessarily for smaller a.
    std::int64_t P3(std::int64_t x, std::int64_t a, int nThreads,
                    const PhiTables& tables);
}
