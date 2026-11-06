#pragma once

#include <cstdint>

namespace PrimeCounting {

    // Sieve counter used by the R driver for 10 <= n < 100000.
    std::int64_t PiPrime(std::int64_t n);

    // Legendre counter for 100000 <= n <= 2^53 - 1; inputs are validated by
    // the R driver. Each call owns its prepared tables and uses PhiParallel;
    // workers share the general Phi evaluator's bounded, thread-safe cache.
    std::int64_t PiLegendre(std::int64_t n, int nThreads = 1,
                                    int maxThreads = 1);

    // Meissel and Lehmer counters for 0 <= n <= 2^53 - 1. Invalid n throws
    // std::out_of_range before table preparation. Each call prepares one
    // immutable square-root table shared by its Phi/P2/P3 phases. Phases run
    // sequentially with at most min(nThreads, maxThreads) requested workers;
    // budgets <= 1 use serial execution. These are internal native entry
    // points; the R method selector and automatic dispatch are separate work.
    std::int64_t PiMeissel(std::int64_t n, int nThreads = 1,
                           int maxThreads = 1);

    std::int64_t PiLehmer(std::int64_t n, int nThreads = 1,
                          int maxThreads = 1);
}
