#pragma once

#include <cstdint>

namespace PrimeCounting {

    // Inclusive sieve counter for 0 <= n <= 2^53 - 1. The R driver uses
    // this for 10 <= n < 100000; larger full-range counts use its existing
    // combinatorial dispatch. No vector of the counted primes is retained.
    std::int64_t PiPrime(std::int64_t n);

    // Count primes in the inclusive interval [minNum, maxNum]. Both bounds
    // must lie in [0, 2^53 - 1]; a reversed interval is empty. Long intervals
    // use adjacent bounded sieve windows. Budgets <= 0 become 1, requests
    // are capped by maxThreads, and small work may use fewer workers.
    std::int64_t PiPrimeRange(std::int64_t minNum, std::int64_t maxNum,
                              int nThreads = 1, int maxThreads = 1,
                              int maxCores = 1);

    // Legendre counter for 1 <= n <= 2^53 - 1; inputs are validated by
    // the R driver. Each call owns its prepared tables and uses PhiParallel;
    // workers share the general Phi evaluator's bounded, thread-safe cache.
    std::int64_t PiLegendre(std::int64_t n, int nThreads = 1,
                                    int maxThreads = 1);

    // Meissel and Lehmer counters for 0 <= n <= 2^53 - 1. Invalid n throws
    // std::out_of_range before table preparation. Each call prepares one
    // immutable square-root table shared by its Phi/P2/P3 phases. Phases run
    // sequentially with at most min(nThreads, maxThreads) requested workers;
    // budgets <= 1 use serial execution. The R driver selects the native
    // counter after validating its public inputs and method argument.
    std::int64_t PiMeissel(std::int64_t n, int nThreads = 1,
                           int maxThreads = 1);

    std::int64_t PiLehmer(std::int64_t n, int nThreads = 1,
                          int maxThreads = 1);
}
