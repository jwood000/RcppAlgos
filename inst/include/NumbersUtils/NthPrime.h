#pragma once

#include <cstdint>

namespace PrimeCounting {

    // pi(2^53 - 1), independently checked with primecount. The limit is on
    // the returned prime, not just the index, so the result remains exact in R.
    constexpr std::int64_t MAX_NTH_PRIME_INDEX = 252252704148404LL;

    // One-based nth prime. Throws std::out_of_range for an unsupported index.
    // Prepared tables are reused across exact counts and the final sieve.
    // nThreads applies to counting; local segmented selection is serial.
    std::int64_t NthPrime(std::int64_t n, int nThreads = 1);
}
