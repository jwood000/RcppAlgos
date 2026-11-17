#pragma once

#include <cstdint>

namespace PrimeCounting {

    class PhiTables;

    enum class PrimeCountMethod {
        Legendre,
        Meissel,
        Lehmer
    };

    // Shared automatic policy. nThreads is the usable budget after the
    // hardware cap; each algorithm can further reduce its worker count.
    // First/warmed-call evidence: .github/development/prime-count-auto/
    inline PrimeCountMethod SelectPrimeCountMethod(std::int64_t n, int nThreads) {
        if (n < 10000000000LL) {
            return PrimeCountMethod::Legendre;
        }

        if (n >= 1000000000000LL && nThreads >= 6) {
            return PrimeCountMethod::Lehmer;
        }

        return PrimeCountMethod::Meissel;
    }

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

    // Legendre counter for 0 <= n <= 2^53 - 1. Each call owns its prepared
    // tables and uses PhiParallel;
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

    // Prepared-table overloads permit repeated counts without rebuilding the
    // prime vector and dense pi table. Require 0 <= n <= 2^53 - 1 and coverage
    // through floor(sqrt(n)); invalid n or coverage throws std::out_of_range.
    // n < 2 returns zero without requiring coverage. Cutoffs depend on n,
    // even when the supplied table covers a larger range. The caller caps
    // nThreads by hardware; budgets <= 1 use serial execution.
    std::int64_t PiLegendre(std::int64_t n, int nThreads,
                            const PhiTables& tables);
    std::int64_t PiMeissel(std::int64_t n, int nThreads,
                           const PhiTables& tables);
    std::int64_t PiLehmer(std::int64_t n, int nThreads,
                          const PhiTables& tables);
    std::int64_t PiAuto(std::int64_t n, int nThreads,
                        const PhiTables& tables);
}
