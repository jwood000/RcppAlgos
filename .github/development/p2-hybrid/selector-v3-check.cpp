// Independent checks of the frozen v3 density refinement. Earlier selector
// fixtures remain separate and unchanged for the comparison experiment.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/Phi.h"
#include "NumbersUtils/ThreadBatch.h"
#include "P2SelectorV3.h"
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
    using P2SelectorV3::Reason;
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

    P2SelectorV3::Decision Choice(const std::vector<std::int64_t>& primes,
                                 std::int64_t x, std::int64_t L, std::int64_t a,
                                 int threads, int candidateMultiplier,
                                 int expectedMultiplier, Reason expectedReason) {
        ++choices;
        const auto range = QueryRange(primes, x, L, a, threads);
        const auto actual = P2SelectorV3::choose(x, L, range.first, range.last,
                                                range.budget, threads, primes);
        Check(actual.multiplier == expectedMultiplier, "unexpected v3 multiplier");
        Check(actual.reason == expectedReason, "unexpected v3 reason");
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

    P2SelectorV3::Decision Kernel(const Oracle& oracle, std::int64_t x,
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
        Check(result == expected, "v3 P2 count mismatch");
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
        const std::vector<std::int64_t> emptyPrimes(1, 0);
        const auto empty = P2SelectorV3::choose(4, 2, 1, 1, 1, 0, emptyPrimes);
        Check(empty.multiplier == 0 && empty.reason == Reason::AtMostOneQuery &&
              empty.limit == 2 && empty.candidateStreamBudget == 0, "empty v3 choice");

        // A very oversized table leaves one or two non-table queries. The
        // two-query span is short but sparse: the unchanged short-full-range
        // rule still accepts it, while one query must retain Phi.
        {
            const std::int64_t x = 100000000;
            const PrimeCounting::PhiTables tables(1000000);
            for (int threads : {-2, 0, 1, 2, 8}) {
                const auto one = Kernel(oracle, x, tables, oracle.pi[97] - 1,
                                         threads, 0, 0, Reason::AtMostOneQuery);
                Check(one.nonTableQueries == 1, "single-query fixture changed");
                const auto two = Kernel(oracle, x, tables, oracle.pi[89] - 1,
                                         threads, 8, 8, Reason::ShortFullRange);
                Check(two.nonTableQueries == 2 && two.candidateSpan == 123595,
                      "short sparse fixture changed");
                Check(two.candidateSpan > 256 * two.candidateStreamQueries,
                      "short-range fixture no longer sparse");
            }
        }

        // Exact density boundaries with full, genuine prime-index ranges.
        // For p=19001, these L values keep the final non-table prime fixed
        // throughout both endpoints of each floor-division plateau:
        //   N=10000 ends at 129967; N=18000 ends at 226669.
        // Neither the oracle nor the expected decisions call the selector.
        struct DensityRange { std::int64_t count, limit; int streamBudget; };
        const DensityRange densityRanges[] = {{10000, 438346, 2}, {18000, 421595, 4}};
        for (const auto& item : densityRanges) {
            const std::int64_t p = 19001;
            const auto a = oracle.pi[p] - 1;
            Check(oracle.primes[a + 1] == p, "density fixture prime");
            const PrimeCounting::PhiTables tables(item.limit);
            for (int delta : {-1, 0, 1}) for (auto remainder : {std::int64_t(0), p - 1}) {
                const auto span = 256 * item.count + delta;
                const auto x = p * (item.limit + span) + remainder;
                Check(x < 100000000000LL, "density fixture left small-input policy");
                for (int threads : {2, 3, 4, 8}) {
                    const bool coverage = threads <= item.streamBudget;
                    const bool accept = coverage && delta <= 0;
                    const Reason reason = !coverage ? Reason::StreamWorkerDeficit :
                        delta > 0 ? Reason::StreamTooSparse : Reason::ParallelCoverage;
                    const auto choice = Kernel(oracle, x, tables, a, threads, 16,
                                               accept ? 16 : 0, reason);
                    Check(choice.candidateStreamQueries == item.count &&
                          choice.nonTableQueries == item.count, "full density query count");
                    Check(choice.candidateSpan == span, "exact full density span");
                    Check(choice.candidateStreamFirst == a + 1, "full density range split");
                }
            }
            // Nonpositive budgets remain serial; sparse candidates do not
            // activate a rule confined to the small-parallel branch.
            const auto x = p * (item.limit + 256 * item.count + 1);
            for (int threads : {-2, 0, 1}) {
                Kernel(oracle, x, tables, a, threads, 8, 8, Reason::Range8);
            }
        }

        // Q=16L is 3/4/5 greater than the last actual target here. At the
        // equality and below-equality cases, substituting Q-L for the actual
        // span would wrongly reject a candidate with sufficient density.
        // An earlier prime query deliberately keeps this range partial.
        {
            const std::int64_t L = 229820, p = 10007, count = 13466;
            const auto firstStream = oracle.pi[p];
            const auto a = firstStream - 2;
            Check(oracle.primes[firstStream - 1] == 9973, "partial density preceding prime");
            const PrimeCounting::PhiTables tables(L);
            for (int delta : {-1, 0, 1}) for (auto remainder : {std::int64_t(0), p - 1}) {
                const auto span = 256 * count + delta;
                const auto x = p * (L + span) + remainder;
                for (int threads : {2, 3, 4, 8}) {
                    const bool coverage = threads <= 3;
                    const bool accept = coverage && delta <= 0;
                    const Reason reason = !coverage ? Reason::StreamWorkerDeficit :
                        delta > 0 ? Reason::StreamTooSparse : Reason::ParallelCoverage;
                    const auto choice = Kernel(oracle, x, tables, a, threads, 16,
                                               accept ? 16 : 0, reason);
                    Check(choice.candidateStreamFirst == firstStream &&
                          choice.candidateStreamQueries == count, "partial density split");
                    Check(choice.candidateSpan == span, "exact partial density span");
                    Check(choice.candidateLimit == 16 * L &&
                          choice.candidateLimit - L > 256 * count,
                          "cap-versus-actual density fixture changed");
                }
            }
        }

        // Preserve the old query guards and their precedence. Wide ranges
        // with enough queries still fail worker coverage before density.
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
                Check(choice.nonTableQueries == n && choice.candidateSpan > 131072,
                      "wide query guard fixture changed");
            }
        }
        for (int L : {47252, 47270}) {
            const PrimeCounting::PhiTables tables(L);
            const auto choice = Kernel(oracle, std::int64_t(L) * L, tables, oracle.pi[300],
                                       1, 8, L == 47252 ? 0 : 8,
                                       L == 47252 ? Reason::TooFewStreamQueries : Reason::Range8);
            Check(choice.candidateStreamQueries == (L == 47252 ? 4095 : 4096),
                  "final stream query guard changed");
        }

        // Rational table sizes are test inputs only: choose receives the
        // actual prepared integer L. These explicit expectations distinguish
        // density rejection from earlier worker rejection as L grows.
        struct TableRatio { int numerator, denominator; Reason two, four, eight; };
        const TableRatio ratios[] = {
            {9, 8, Reason::ParallelCoverage, Reason::ParallelCoverage, Reason::StreamWorkerDeficit},
            {6, 5, Reason::StreamTooSparse, Reason::StreamTooSparse, Reason::StreamWorkerDeficit},
            {5, 4, Reason::StreamTooSparse, Reason::StreamTooSparse, Reason::StreamWorkerDeficit},
            {3, 2, Reason::StreamTooSparse, Reason::StreamTooSparse, Reason::StreamWorkerDeficit},
            {2, 1, Reason::StreamTooSparse, Reason::StreamTooSparse, Reason::StreamTooSparse},
            {3, 1, Reason::StreamTooSparse, Reason::StreamTooSparse, Reason::StreamTooSparse}
        };
        for (const auto& item : ratios) {
            const std::int64_t x = 70000000000LL;
            const auto L = std::int64_t(Root(x)) * item.numerator / item.denominator;
            const auto a = oracle.pi[10000];
            const PrimeCounting::PhiTables tables(L);
            for (int threads : {-2, 0, 1}) {
                Kernel(oracle, x, tables, a, threads, 8, 8, Reason::Range8);
            }
            for (int threads : {2, 4, 8}) {
                const Reason reason = threads == 2 ? item.two : threads == 4 ? item.four : item.eight;
                const bool accept = reason == Reason::ParallelCoverage;
                Kernel(oracle, x, tables, a, threads, 16, accept ? 16 : 0, reason);
            }
        }
        {
            const std::int64_t x = 70000000000LL;
            const auto L = 2 * Root(x);
            const PrimeCounting::PhiTables tables(L);
            const auto full = QueryRange(oracle.primes, x, L, 0, 4);
            for (int n : {8191, 8192, 8193}) {
                const auto a = full.last - 1 - n;
                Kernel(oracle, x, tables, a, 4, 16, 0, Reason::StreamWorkerDeficit);
                Check(P2HybridMetrics::current().phiWorkers == (n < 8192 ? 3 : 4),
                      "original effective-worker cliff changed");
            }
        }

        // Density is confined to x<1e11 with a parallel Phi budget. The same
        // sparse geometry on/above that boundary must keep v2's large-input
        // choice, including cases that do not satisfy small-input parity.
        for (std::int64_t x : {99999999999LL, 100000000000LL, 100000000001LL}) {
            const auto L = 2 * Root(x);
            const PrimeCounting::PhiTables tables(L);
            const auto a = oracle.pi[15000];
            for (int threads : {1, 2, 4}) {
                const bool small = x < 100000000000LL;
                const int k = small && threads == 1 ? 8 : 16;
                const bool reject = small && threads > 1;
                const Reason reason = reject ? Reason::StreamTooSparse :
                    small ? Reason::Range8 : Reason::Range16;
                const auto choice = Kernel(oracle, x, tables, a, threads, k, reject ? 0 : k, reason);
                Check(choice.candidateSpan > 256 * choice.candidateStreamQueries,
                      "large-input sparse exception fixture changed");
            }
        }

        // Direct large-input choices retain all earlier caps/promotions. The
        // common checker separately exercises complete Meissel/Lehmer counts.
        struct Cap { std::int64_t x; int serial, parallel; Reason parallelReason; };
        const Cap caps[] = {
            {999999999999LL, 16, 16, Reason::Range16},
            {1000000000000LL, 16, 32, Reason::ParallelPromotion},
            {2999999999999LL, 16, 32, Reason::ParallelPromotion},
            {3000000000000LL, 32, 32, Reason::Range32},
            {29999999999999LL, 32, 32, Reason::Range32},
            {30000000000000LL, 64, 64, Reason::Range64}
        };
        for (const auto& item : caps) {
            const Reason serialReason = item.serial == 16 ? Reason::Range16 :
                item.serial == 32 ? Reason::Range32 : Reason::Range64;
            Choice(oracle.primes, item.x, Root(item.x), 0, 1,
                   item.serial, item.serial, serialReason);
            Choice(oracle.primes, item.x, Root(item.x), 0, 8,
                   item.parallel, item.parallel, item.parallelReason);
        }

        // Genuine complete prime coverage at MAX_X, without a dense pi table
        // or a giant P2 run. Wide stream windows are checked separately.
        {
            constexpr std::int64_t maximum = 9007199254740991LL;
            const auto L = Root(maximum);
            const auto primes = PrimeList(L);
            Check(L == 94906265 && primes.size() > 5000000, "maximum prime coverage");
            const auto cutoff = maximum / (40 * std::int64_t(L));
            const auto a = std::upper_bound(primes.begin(), primes.end(), cutoff) - primes.begin() - 1;
            const auto uncapped = Choice(primes, maximum, L, 0, 8, 64, 64, Reason::Range64);
            Check(uncapped.limit == 64 * std::int64_t(L), "maximum uncapped endpoint");
            Check(uncapped.candidateSpan > 256 * uncapped.candidateStreamQueries,
                  "MAX_X density exception fixture changed");
            for (int threads : {-2, 0, 8}) {
                const auto capped = Choice(primes, maximum - (threads == 0), L, a,
                                           threads, 64, 64, Reason::Range64);
                Check(capped.limit < 64 * std::int64_t(L) &&
                      capped.limit > std::numeric_limits<std::int32_t>::max(),
                      "maximum capped 64-bit endpoint");
                Check(capped.streamFirst == a + 1, "maximum capped full range");
            }
        }

        // Concurrent callers share tables/cache while selecting opposite
        // sides of the exact density boundary. Repeating after join also
        // verifies the same decisions/results with warmed shared state.
        {
            const std::int64_t L = 438346, p = 19001, count = 10000;
            const auto a = oracle.pi[p] - 1;
            const PrimeCounting::PhiTables tables(L);
            const std::int64_t xs[2] = {p * (L + 256 * count), p * (L + 256 * count + 1)};
            const std::int64_t expected[2] = {oracle.pairs(xs[0], a), oracle.pairs(xs[1], a)};
            for (int repetition = 0; repetition < 2; ++repetition) {
                std::int64_t sums[2] = {0, 0}, streams[2] = {0, 0};
                PrimeCounting::ThreadBatch callers(2);
                for (int slot = 0; slot < 2; ++slot) callers.launch([&, slot] {
                    sums[slot] = PrimeCounting::P2(xs[slot], a, 2, tables);
                    streams[slot] = P2HybridMetrics::current().streamQueries;
                });
                callers.finish();
                Check(sums[0] == expected[0] && sums[1] == expected[1], "concurrent v3 counts");
                Check(streams[0] == count && streams[1] == 0, "concurrent v3 density choices");
            }
        }

        std::cout << "PASS selector_v3_checks=" << checks << " choice_cases=" << choices
                  << " kernel_cases=" << kernels << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
