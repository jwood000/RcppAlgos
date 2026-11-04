#pragma once

#include <cstdint>

namespace PrimeCounting {

    // Sieve counter used by the R driver for 10 <= n < 100000.
    std::int64_t PiPrime(std::int64_t n);

    // Legendre counter for 100000 <= n <= 2^53 - 1; inputs are validated by
    // the R driver. Each call owns its prepared tables; workers share the
    // general Phi evaluator's bounded, thread-safe cache.
    std::int64_t PiLegendre(std::int64_t n, int nThreads = 1,
                                    int maxThreads = 1);
}
