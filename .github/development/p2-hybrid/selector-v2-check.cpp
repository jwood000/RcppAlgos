// Independent checks of the frozen, interval-aware v2 experiment. Keep the
// v1 fixture separate so evaluating this policy cannot alter its baseline.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/Phi.h"
#include "NumbersUtils/ThreadBatch.h"
#include "P2SelectorV2.h"
#include "hybrid_metrics.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace PrimeCounting {
    std::int64_t P2Reference(std::int64_t, std::int64_t, int, const PhiTables&);
}

namespace {
    using P2SelectorV2::Reason;
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
            const auto root = Root(x);
            Check(root < static_cast<int>(pi.size()), "short independent root coverage");
            const auto b = pi[root];
            if (a < b) {
                Check(x / primes[a + 1] < static_cast<std::int64_t>(pi.size()),
                      "short independent quotient coverage");
            }
            std::int64_t result = 0;
            for (auto i = a + 1; i <= b; ++i) {
                result += pi[x / primes[i]] - i + 1;
            }
            return result;
        }
    };

    // MAX_X choice checks need genuine primes through its integer root, but
    // neither a dense pi table nor a complete P2 calculation. This separate
    // odd-only sieve uses about 6 MB of bits plus the retained prime vector.
    std::vector<std::int64_t> PrimeList(int limit) {
        std::vector<bool> composite(static_cast<std::size_t>(limit / 2 + 1), false);
        for (int p = 3; p <= limit / p; p += 2) if (!composite[p / 2]) {
            for (int n = p * p; n <= limit; n += 2 * p) composite[n / 2] = true;
        }
        std::vector<std::int64_t> primes = {0, 2};
        for (int n = 3; n <= limit; n += 2) if (!composite[n / 2]) primes.push_back(n);
        return primes;
    }

    struct Range {
        std::int64_t first, last, b;
        std::size_t budget;
    };

    Range QueryRange(const std::vector<std::int64_t>& primes, std::int64_t x,
                     std::int64_t L, std::int64_t a, int threads) {
        Check(Root(x) <= L, "invalid prepared limit in fixture");
        const auto b = std::upper_bound(primes.begin(), primes.end(), Root(x)) -
            primes.begin() - 1;
        Range range = {a + 1, b + 1, b, 1};
        // Direct integer quotients independently establish the range; do not
        // reuse the kernel's floor(x/(Q+1)) binary-search transformation.
        while (range.last > range.first && x / primes[range.last - 1] <= L) --range.last;
        range.budget = static_cast<std::size_t>(std::min<std::int64_t>(
            std::max(1, threads), std::max<std::int64_t>(1, (range.last - range.first) / 2048)));
        return range;
    }

    P2SelectorV2::Decision Choice(const std::vector<std::int64_t>& primes,
                                 std::int64_t x, std::int64_t L, std::int64_t a,
                                 int threads, int candidateMultiplier,
                                 int expectedMultiplier, Reason expectedReason) {
        ++choices;
        const auto range = QueryRange(primes, x, L, a, threads);
        const auto actual = P2SelectorV2::choose(x, L, range.first, range.last,
                                                range.budget, threads, primes);
        Check(actual.multiplier == expectedMultiplier, "unexpected v2 multiplier");
        Check(actual.reason == expectedReason, "unexpected v2 reason");
        // These fixtures use L <= 94906265 and k <= 64, so independently
        // multiplying before capping is exact and well within int64_t.
        const auto candidateLimit = candidateMultiplier ?
            std::min<std::int64_t>(L * candidateMultiplier, x / primes[range.first]) : L;
        auto candidateFirst = range.first;
        while (candidateFirst < range.last && x / primes[candidateFirst] > candidateLimit) {
            ++candidateFirst;
        }
        const auto candidateQueries = range.last - candidateFirst;
        const auto span = candidateQueries ? x / primes[candidateFirst] - L : 0;
        const auto streamBudget = candidateQueries ? std::min<std::int64_t>(
            std::max(1, threads), std::min(candidateQueries, std::max<std::int64_t>(1, span / 1048576))) : 0;
        Check(actual.candidateLimit == candidateLimit, "candidate endpoint mismatch");
        Check(actual.candidateStreamFirst == candidateFirst, "candidate split mismatch");
        Check(actual.candidateStreamQueries == candidateQueries, "candidate query count mismatch");
        Check(actual.candidateSpan == span, "candidate actual span mismatch");
        Check(actual.candidateStreamBudget == static_cast<std::size_t>(streamBudget),
              "candidate stream worker budget mismatch");
        const auto expectedLimit = expectedMultiplier ? candidateLimit : L;
        const auto expectedFirst = expectedMultiplier ? candidateFirst : range.last;
        Check(actual.limit == expectedLimit, "selected endpoint mismatch");
        Check(actual.streamFirst == expectedFirst, "selected split mismatch");
        Check(actual.nonTableQueries == range.last - range.first, "non-table query count mismatch");
        Check(actual.phiQueries == expectedFirst - range.first, "Phi query count mismatch");
        Check(actual.streamQueries == range.last - expectedFirst, "stream query count mismatch");
        Check(actual.phiQueries + actual.streamQueries == actual.nonTableQueries,
              "selector partition does not cover the non-table range");
        if (candidateQueries) {
            Check(Root(x / primes[candidateFirst]) <= L, "stream base-prime coverage");
            Check(x / primes[candidateFirst] <= candidateLimit, "first stream target exceeds cap");
            if (candidateFirst > range.first) {
                Check(x / primes[candidateFirst - 1] > candidateLimit,
                      "Phi prefix overlaps the stream interval");
            }
        }
        return actual;
    }

    P2SelectorV2::Decision Kernel(const Oracle& oracle, std::int64_t x,
                                 const PrimeCounting::PhiTables& tables,
                                 std::int64_t a, int threads, int candidateMultiplier,
                                 int multiplier, Reason reason, bool independent = true) {
        ++kernels;
        const auto choice = Choice(oracle.primes, x, tables.limit(), a, threads,
                                   candidateMultiplier, multiplier, reason);
        const auto range = QueryRange(oracle.primes, x, tables.limit(), a, threads);
        const auto expected = independent ? oracle.pairs(x, a) :
            PrimeCounting::P2Reference(x, a, 1, tables);
        const auto result = PrimeCounting::P2(x, a, threads, tables);
        const auto& metrics = P2HybridMetrics::current();
        Check(result == expected, "v2 P2 count mismatch");
        Check(metrics.streamLimit == choice.limit, "kernel selected endpoint mismatch");
        Check(metrics.phiQueries == choice.phiQueries, "kernel Phi partition mismatch");
        Check(metrics.streamQueries == choice.streamQueries, "kernel stream partition mismatch");
        Check(metrics.directQueries == range.b + 1 - range.last, "kernel direct partition mismatch");
        const auto phiWorkers = choice.phiQueries ?
            std::min<std::size_t>(range.budget, static_cast<std::size_t>(choice.phiQueries)) : 0;
        Check(metrics.phiWorkers == phiWorkers, "original Phi worker budget changed");
        const auto streamBudget = choice.streamQueries ? choice.candidateStreamBudget : 0;
        Check(metrics.intervals.size() == streamBudget, "kernel stream interval budget mismatch");
        Check(metrics.streamWorkers <= streamBudget, "kernel stream worker budget exceeded");
        Check((metrics.streamWorkers == 0) == (choice.streamQueries == 0),
              "kernel stream workers do not match selected work");
        return choice;
    }
}

int main() {
    try {
        const Oracle oracle;
        // No prime access is permitted for an empty range. One query keeps
        // Phi even when streaming would otherwise qualify as a short range.
        const std::vector<std::int64_t> emptyPrimes(1, 0);
        const auto empty = P2SelectorV2::choose(4, 2, 1, 1, 1, 0, emptyPrimes);
        Check(empty.multiplier == 0 && empty.reason == Reason::AtMostOneQuery &&
              empty.limit == 2 && empty.candidateStreamBudget == 0, "empty v2 choice");
        {
            const std::int64_t x = 1000000;
            const PrimeCounting::PhiTables tables(1000);
            for (int threads : {-2, 0, 1, 2, 8}) {
                Kernel(oracle, x, tables, oracle.pi[1000] - 1, threads, 0, 0,
                       Reason::AtMostOneQuery);
                const auto choice = Kernel(oracle, x, tables, oracle.pi[1000] - 2,
                                           threads, 8, 8, Reason::ShortFullRange);
                Check(choice.nonTableQueries == 2, "two-query fixture changed");
            }
        }

        // Exact span and floor-division boundaries, with genuine prime-index
        // ranges. The 1999 cutoff leaves one Phi query; 2003 covers the whole
        // range. Q=8L=149800 is slightly larger than the actual last target.
        {
            const std::int64_t L = 18725, p = 2003;
            Check(oracle.primes[oracle.pi[1999]] == 1999 &&
                  oracle.primes[oracle.pi[1999] + 1] == p, "adjacent-prime fixture changed");
            const PrimeCounting::PhiTables tables(L);
            for (int span : {131071, 131072, 131073}) {
                for (int remainder : {0, 2002}) {
                    const auto x = p * (L + span) + remainder;
                    for (bool full : {false, true}) {
                        const auto a = oracle.pi[1999] - (full ? 0 : 1);
                        for (int threads : {-2, 0, 1, 2, 8}) {
                            const bool accept = full && span <= 131072;
                            const auto choice = Kernel(oracle, x, tables, a, threads, 8,
                                accept ? 8 : 0, accept ? Reason::ShortFullRange : Reason::TooFewQueries);
                            Check(choice.candidateSpan == span, "exact short-span fixture changed");
                            Check(choice.nonTableQueries < 4096, "short-span guard fixture changed");
                            Check((choice.candidateStreamFirst == a + 1) == full,
                                  "full/partial range fixture changed");
                        }
                    }
                }
            }
        }

        // A short full range bypasses the old guard and the small-parallel
        // worker comparison, including exactly 4095/4096/4097 queries.
        {
            const std::int64_t x = 10000000000LL;
            const PrimeCounting::PhiTables tables(100000);
            for (int n : {4095, 4096, 4097}) for (int threads : {-2, 0, 1, 2, 8}) {
                const int k = n >= 4096 && threads > 1 ? 16 : 8;
                const auto choice = Kernel(oracle, x, tables, oracle.pi[100000] - n,
                                           threads, k, k, Reason::ShortFullRange);
                Check(choice.nonTableQueries == n, "short query-count fixture changed");
            }
        }

        // Wide ranges still obey both 4096-query guards. Oversized tables
        // change the actual non-table endpoint independently of pi(sqrt(x)).
        {
            const std::int64_t x = 10000000000LL;
            const PrimeCounting::PhiTables tables(200000);
            const auto full = QueryRange(oracle.primes, x, tables.limit(), 0, 1);
            for (int n : {4095, 4096, 4097}) for (int threads : {1, 2, 8}) {
                const auto a = full.last - 1 - n;
                const bool parallel = n >= 4096 && threads > 1;
                const int k = parallel ? 16 : 8;
                const Reason reason = n < 4096 ? Reason::TooFewQueries :
                    parallel ? Reason::StreamWorkerDeficit : Reason::Range8;
                const auto choice = Kernel(oracle, x, tables, a, threads, k,
                                           n >= 4096 && !parallel ? 8 : 0, reason);
                Check(choice.candidateSpan > 131072 && choice.nonTableQueries == n,
                      "wide query-count fixture changed");
            }
        }
        for (int L : {47252, 47270}) {
            const PrimeCounting::PhiTables tables(L);
            const auto choice = Kernel(oracle, std::int64_t(L) * L, tables, oracle.pi[300],
                                       1, 8, L == 47252 ? 0 : 8,
                                       L == 47252 ? Reason::TooFewStreamQueries : Reason::Range8);
            Check(choice.candidateStreamQueries == (L == 47252 ? 4095 : 4096),
                  "final stream guard fixture changed");
            Check(choice.candidateSpan > 131072, "wide final stream guard fixture changed");
        }

        // 2- and 4-worker numerical-span boundaries. Test both ends of each
        // floor quotient's plateau, actual worker parity, and a larger Phi
        // budget that must reject the same candidate. These are real ranges
        // with all quotients within the independent 9e6 oracle.
        for (int workerCliff : {2, 4}) {
            const std::int64_t L = workerCliff == 2 ? 150000 : 300000;
            const auto p = oracle.primes[oracle.pi[workerCliff == 2 ? 10000 : 19000] + 1];
            const auto a = oracle.pi[p] - 1;
            const PrimeCounting::PhiTables tables(L);
            for (int delta : {-1, 0, 1}) for (auto remainder : {std::int64_t(0), p - 1}) {
                const auto span = std::int64_t(workerCliff) * 1048576 + delta;
                const auto x = p * (L + span) + remainder;
                Check(x < 100000000000LL, "worker-cliff fixture left the small-input policy");
                for (int threads : {2, 4, 8}) {
                    const bool accept = threads <= workerCliff &&
                        (threads < workerCliff || delta >= 0);
                    const auto choice = Kernel(oracle, x, tables, a, threads, 16,
                        accept ? 16 : 0, accept ? Reason::ParallelCoverage : Reason::StreamWorkerDeficit);
                    Check(choice.candidateSpan == span, "exact worker span fixture changed");
                }
            }
            const auto x = p * (L + std::int64_t(workerCliff) * 1048576);
            for (int threads : {-2, 0, 1}) {
                Kernel(oracle, x, tables, a, threads, 8, 8, Reason::Range8);
            }
        }

        // Preserve v1's large-input caps and effective-worker promotion at
        // each boundary. Full-range choices expose uncapped Q; kernel calls
        // use smaller, independently countable quotient ranges.
        struct Policy { std::int64_t x; int serial, parallel; Reason serialReason, parallelReason; };
        const Policy policy[] = {
            {99999999999LL, 8, 0, Reason::Range8, Reason::StreamWorkerDeficit},
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
            const int parallelCandidate = item.parallel ? item.parallel : 16;
            Choice(oracle.primes, item.x, L, 0, 1, item.serial, item.serial, item.serialReason);
            Choice(oracle.primes, item.x, L, 0, 8, parallelCandidate, item.parallel, item.parallelReason);
            const PrimeCounting::PhiTables tables(L);
            const auto a = oracle.pi[2 * L / 3];
            Kernel(oracle, item.x, tables, a, 1, item.serial, item.serial, item.serialReason);
            Kernel(oracle, item.x, tables, a, 8, parallelCandidate, item.parallel, item.parallelReason);
        }
        {
            const std::int64_t x = 1000000000000LL;
            const PrimeCounting::PhiTables tables(1000000);
            for (int n : {8191, 8192, 8193}) for (int threads : {3, 4, 8}) {
                const int k = n >= 8192 && threads >= 4 ? 32 : 16;
                Kernel(oracle, x, tables, oracle.pi[1000000] - n, threads,
                       k, k, Reason::ShortFullRange);
            }
            // Two explicitly differential checks expose the changed Q in
            // actual kernel work, rather than just the reported multiplier.
            Kernel(oracle, x, tables, oracle.pi[10000], 3, 16, 16, Reason::Range16, false);
            Kernel(oracle, x, tables, oracle.pi[10000], 4, 32, 32, Reason::ParallelPromotion, false);
        }

        // Exercise the supported maximum without allocating a full pi table
        // or attempting to stream billions of numbers. The complete genuine
        // prime vector keeps all chooser preconditions valid. The existing
        // stream-window fixture independently checks these wide coordinates.
        {
            constexpr std::int64_t maximum = 9007199254740991LL;
            const auto L = Root(maximum);
            const auto primes = PrimeList(L);
            Check(L == 94906265, "maximum integer root");
            Check(primes.back() <= L && primes.size() > 5000000, "maximum prime-list coverage");
            const auto cutoff = maximum / (40 * std::int64_t(L));
            const auto a = std::upper_bound(primes.begin(), primes.end(), cutoff) - primes.begin() - 1;
            const auto uncapped = Choice(primes, maximum, L, 0, 8, 64, 64, Reason::Range64);
            Check(uncapped.limit == 64 * std::int64_t(L), "maximum uncapped endpoint");
            for (int threads : {-2, 0, 8}) {
                const auto x = maximum - (threads == 0);
                const auto capped = Choice(primes, x, L, a, threads, 64, 64, Reason::Range64);
                Check(capped.limit < 64 * std::int64_t(L) &&
                      capped.limit > std::numeric_limits<std::int32_t>::max(),
                      "maximum capped 64-bit endpoint");
                Check(capped.streamFirst == a + 1, "maximum capped full range");
            }
        }

        // Shared immutable context, adjacent spans, opposite policy choices.
        // Per-caller metrics must stay separate while worker groups overlap.
        {
            const std::int64_t L = 300000;
            const auto p = oracle.primes[oracle.pi[19000] + 1];
            const auto a = oracle.pi[p] - 1;
            const PrimeCounting::PhiTables tables(L);
            std::int64_t sums[2] = {0, 0}, streams[2] = {0, 0};
            std::int64_t xs[2] = {p * (L + 4 * 1048576LL - 1), p * (L + 4 * 1048576LL)};
            const std::int64_t expected[2] = {oracle.pairs(xs[0], a), oracle.pairs(xs[1], a)};
            PrimeCounting::ThreadBatch callers(2);
            for (int slot = 0; slot < 2; ++slot) callers.launch([&, slot] {
                sums[slot] = PrimeCounting::P2(xs[slot], a, 4, tables);
                streams[slot] = P2HybridMetrics::current().streamQueries;
            });
            callers.finish();
            Check(sums[0] == expected[0] && sums[1] == expected[1], "concurrent v2 P2 results");
            Check(streams[0] == 0 && streams[1] >= 4096, "concurrent opposite v2 choices");
        }

        std::cout << "PASS selector_v2_checks=" << checks << " choice_cases=" << choices
                  << " kernel_cases=" << kernels << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
