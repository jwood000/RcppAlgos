#pragma once

#include <cstdint>

namespace PrimeCounting {

    // pi(2^53 - 1), independently checked with primecount. The limit is on
    // the returned prime, not just the index, so the result remains exact in R.
    constexpr std::int64_t MAX_NTH_PRIME_INDEX = 252252704148404LL;

    // One-based nth prime. Throws std::out_of_range for an unsupported index.
    // Prepared tables are reused across at most two whole-prefix counts,
    // then released before interval counting and bounded native selection.
    // The two-argument entry uses nThreads as its usable hardware budget.
    std::int64_t NthPrime(std::int64_t n, int nThreads = 1);

    // Full hardware budgets preserve the same sieve policy as primeSieve.
    // Requests are capped by maxThreads; nonpositive budgets become one.
    std::int64_t NthPrime(std::int64_t n, int nThreads,
                          int maxThreads, int maxCores);
}
