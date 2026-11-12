// Frozen-policy checks and independent integration checks for the temporary
// P2 selector. Compile the selected kernel with P2_HYBRID_SELECTOR and
// P2_HYBRID_MEMORY; timings must use a different, uninstrumented binary.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/Phi.h"
#include "NumbersUtils/ThreadBatch.h"
#include "P2Selector.h"
#include "hybrid_metrics.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace PrimeCounting {
    std::int64_t P2Reference(std::int64_t, std::int64_t, int, const PhiTables&);
}

namespace {
    std::size_t checks = 0, choices = 0, kernels = 0;

    void Check(bool condition, const char* label) {
        ++checks;
        if (!condition) throw std::runtime_error(label);
    }

    int Root(std::int64_t x) {
        int low = 0, high = static_cast<int>(std::min<std::int64_t>(x, 94906265));
        while (low < high) {
            const int middle = low + (high - low + 1) / 2;
            if (middle <= x / middle) low = middle;
            else high = middle - 1;
        }
        return low;
    }

    struct Oracle {
        std::vector<int> pi;
        std::vector<std::int64_t> primes;

        Oracle() : pi(9000001, 0), primes(1, 0) {
            std::vector<unsigned char> sieve(pi.size(), 1);
            sieve[0] = sieve[1] = 0;
            const int maximum = static_cast<int>(pi.size()) - 1;
            for (int p = 2; p <= maximum / p; ++p) if (sieve[p]) {
                for (int n = p * p; n <= maximum; n += p) sieve[n] = 0;
            }
            for (int n = 1; n <= maximum; ++n) {
                pi[n] = pi[n - 1] + sieve[n];
                if (sieve[n]) primes.push_back(n);
            }
        }

        std::int64_t pairs(std::int64_t x, std::int64_t a) const {
            std::int64_t result = 0;
            for (auto i = a + 1; primes[i] <= x / primes[i]; ++i) {
                const auto q = x / primes[i];
                Check(q < static_cast<std::int64_t>(pi.size()), "short independent pi oracle");
                result += pi[q] - i + 1;
            }
            return result;
        }
    };

    struct Range {
        std::int64_t first, last, b;
        std::size_t budget;
    };

    Range QueryRange(const Oracle& oracle, std::int64_t x, std::int64_t L,
                     std::int64_t a, int threads) {
        const std::int64_t b = oracle.pi[Root(x)];
        Range range = {a + 1, b + 1, b, 1};
        // Linear endpoint checks are independent of the kernel's binary
        // search and floor(x/(L+1)) transformation.
        while (range.last > range.first && x / oracle.primes[range.last - 1] <= L) {
            --range.last;
        }
        range.budget = static_cast<std::size_t>(std::min<std::int64_t>(
            std::max(1, threads), std::max<std::int64_t>(1, (range.last - range.first) / 2048)));
        return range;
    }

    P2Selector::Decision Choice(const Oracle& oracle, std::int64_t x,
                                std::int64_t L, std::int64_t a, int threads,
                                int expectedMultiplier, P2Selector::Reason expectedReason) {
        ++choices;
        const auto range = QueryRange(oracle, x, L, a, threads);
        const auto actual = P2Selector::choose(x, L, range.first, range.last,
                                              range.budget, oracle.primes);
        Check(actual.multiplier == expectedMultiplier, "unexpected frozen multiplier");
        Check(actual.reason == expectedReason, "unexpected frozen reason");
        const auto expectedLimit = expectedMultiplier ?
            std::min<std::int64_t>(L * expectedMultiplier, x / oracle.primes[range.first]) : L;
        auto expectedFirst = range.first;
        while (expectedFirst < range.last && x / oracle.primes[expectedFirst] > expectedLimit) {
            ++expectedFirst;
        }
        Check(actual.limit == expectedLimit, "unexpected endpoint");
        Check(actual.streamFirst == expectedFirst, "unexpected first streamed index");
        Check(actual.nonTableQueries == range.last - range.first, "non-table query count");
        Check(actual.phiQueries == expectedFirst - range.first, "Phi query count");
        Check(actual.streamQueries == range.last - expectedFirst, "stream query count");
        Check(actual.phiQueries + actual.streamQueries == actual.nonTableQueries,
              "selector partition does not cover non-table queries");
        return actual;
    }

    void Kernel(const Oracle& oracle, std::int64_t x, const PrimeCounting::PhiTables& tables,
                std::int64_t a, int threads, int multiplier, P2Selector::Reason reason,
                bool independent = true) {
        ++kernels;
        const auto choice = Choice(oracle, x, tables.limit(), a, threads, multiplier, reason);
        const auto range = QueryRange(oracle, x, tables.limit(), a, threads);
        const auto expected = independent ? oracle.pairs(x, a) :
            PrimeCounting::P2Reference(x, a, 1, tables);
        const auto result = PrimeCounting::P2(x, a, threads, tables);
        const auto& metrics = P2HybridMetrics::current();
        Check(result == expected, "selected P2 count mismatch");
        Check(metrics.streamLimit == choice.limit, "kernel did not use selected endpoint");
        Check(metrics.phiQueries == choice.phiQueries, "kernel Phi partition mismatch");
        Check(metrics.streamQueries == choice.streamQueries, "kernel stream partition mismatch");
        Check(metrics.directQueries == range.b + 1 - range.last, "kernel direct partition mismatch");
        const auto expectedPhiWorkers = choice.phiQueries ?
            std::min<std::size_t>(range.budget, static_cast<std::size_t>(choice.phiQueries)) : 0;
        Check(metrics.phiWorkers == expectedPhiWorkers, "original Phi budget was not retained");
        Check(metrics.streamWorkers <= static_cast<std::uint64_t>(std::max(1, threads)),
              "stream worker budget exceeded");
        Check((metrics.streamWorkers == 0) == (choice.streamQueries == 0),
              "stream workers do not match selected work");
    }
}

int main() {
    try {
        using P2Selector::Reason;
        const Oracle oracle;

        // Empty/short ranges must return before indexing primes[first].
        const std::vector<std::int64_t> emptyPrimes(1, 0);
        const auto empty = P2Selector::choose(4, 2, 1, 1, 1, emptyPrimes);
        Check(empty.multiplier == 0 && empty.reason == Reason::TooFewQueries &&
              empty.limit == 2 && empty.streamQueries == 0, "empty selector range");

        // Uncapped ranges expose the multiplier itself at every x boundary.
        struct Policy { std::int64_t x; int serial, parallel; Reason serialReason, parallelReason; };
        const Policy policy[] = {
            {99999999999LL, 8, 0, Reason::Range8, Reason::SmallParallel},
            {100000000000LL, 16, 16, Reason::Range16, Reason::Range16},
            {100000000001LL, 16, 16, Reason::Range16, Reason::Range16},
            {999999999999LL, 16, 16, Reason::Range16, Reason::Range16},
            {1000000000000LL, 16, 32, Reason::Range16, Reason::ParallelPromotion},
            {1000000000001LL, 16, 32, Reason::Range16, Reason::ParallelPromotion},
            {2999999999999LL, 16, 32, Reason::Range16, Reason::ParallelPromotion},
            {3000000000000LL, 32, 32, Reason::Range32, Reason::Range32},
            {3000000000001LL, 32, 32, Reason::Range32, Reason::Range32},
            {29999999999999LL, 32, 32, Reason::Range32, Reason::Range32},
            {30000000000000LL, 64, 64, Reason::Range64, Reason::Range64},
            {30000000000001LL, 64, 64, Reason::Range64, Reason::Range64}
        };
        for (const auto& item : policy) {
            const auto L = Root(item.x);
            Choice(oracle, item.x, L, 0, 1, item.serial, item.serialReason);
            Choice(oracle, item.x, L, 0, 8, item.parallel, item.parallelReason);
            // Close-to-root cutoffs keep every quotient below9e6, providing
            // independent result checks on both sides without huge sweeps.
            const auto a = oracle.pi[L] - 8192;
            const PrimeCounting::PhiTables tables(L);
            Kernel(oracle, item.x, tables, a, 1, item.serial, item.serialReason);
            Kernel(oracle, item.x, tables, a, 8, item.parallel, item.parallelReason);
        }

        // Non-table query guards and requested-vs-effective worker budgets.
        {
            const std::int64_t x = 10000000000LL;
            const int L = Root(x), b = oracle.pi[L];
            const PrimeCounting::PhiTables tables(L);
            for (int n : {4095, 4096, 4097}) {
                for (int threads : {-2, 0, 1, 2, 8}) {
                    const bool tooFew = n < 4096;
                    const bool parallel = threads > 1 && !tooFew;
                    Kernel(oracle, x, tables, b - n, threads,
                           tooFew || parallel ? 0 : 8,
                           tooFew ? Reason::TooFewQueries :
                           parallel ? Reason::SmallParallel : Reason::Range8);
                }
            }
        }
        {
            const std::int64_t x = 1000000000000LL;
            const int L = Root(x), b = oracle.pi[L];
            const PrimeCounting::PhiTables tables(L);
            for (int n : {8191, 8192, 8193}) {
                for (int threads : {3, 4, 8}) {
                    const bool promote = n >= 8192 && threads >= 4;
                    Kernel(oracle, x, tables, b - n, threads, promote ? 32 : 16,
                           promote ? Reason::ParallelPromotion : Reason::Range16);
                }
            }
            // Here Q is not capped: the kernel's use of3vs4 effective Phi
            // workers must also change its actual streamed endpoint.
            const auto a = oracle.pi[10000];
            Kernel(oracle, x, tables, a, 3, 16, Reason::Range16, false);
            Kernel(oracle, x, tables, a, 4, 32, Reason::ParallelPromotion, false);
        }

        // Exact final-stream guard: these endpoints contain4095/4096 eligible
        // streamed queries under8L. Cutoff pi(300) keeps q inside the oracle.
        for (int L : {47252, 47270}) {
            const std::int64_t x = std::int64_t(L) * L;
            const PrimeCounting::PhiTables tables(L);
            const auto a = oracle.pi[300];
            Kernel(oracle, x, tables, a, 1, L == 47252 ? 0 : 8,
                   L == 47252 ? Reason::TooFewStreamQueries : Reason::Range8);
        }

        // Oversized tables move the direct/non-table boundary independently of
        // a; the selector must use the actual remaining range and cap its Q.
        {
            const std::int64_t x = 10000000000LL;
            const int L = 2 * Root(x);
            const PrimeCounting::PhiTables tables(L);
            const auto full = QueryRange(oracle, x, L, 0, 1);
            for (int n : {4095, 4096}) {
                const auto a = full.last - 1 - n;
                Kernel(oracle, x, tables, a, 1, n == 4095 ? 0 : 8,
                       n == 4095 ? Reason::TooFewQueries : Reason::Range8);
            }
        }

        // Opposite decisions can coexist while sharing the immutable context
        // and Phi cache. Diagnostics are calling-thread-local, including their
        // worker-written interval slots; collect them before callers finish.
        {
            const std::int64_t x = 10000000000LL;
            const int L = Root(x);
            const PrimeCounting::PhiTables tables(L);
            const auto a = oracle.pi[L] - 8192;
            const auto expected = oracle.pairs(x, a);
            std::int64_t sums[2] = {0, 0}, streams[2] = {0, 0};
            PrimeCounting::ThreadBatch callers(2);
            for (int slot = 0; slot < 2; ++slot) {
                callers.launch([&, slot] {
                    sums[slot] = PrimeCounting::P2(x, a, slot ? 8 : 1, tables);
                    streams[slot] = P2HybridMetrics::current().streamQueries;
                });
            }
            callers.finish();
            Check(sums[0] == expected && sums[1] == expected, "concurrent selected counts");
            Check(streams[0] >= 4096 && streams[1] == 0, "concurrent selector decisions");
        }

        std::cout << "PASS selector_checks=" << checks << " choice_cases=" << choices
                  << " kernel_cases=" << kernels << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
