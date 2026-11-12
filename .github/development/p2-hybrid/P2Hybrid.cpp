// Experimental drop-in P2 implementation; never compiled by the R package.
// Build variants with -DP2_HYBRID_MULTIPLIER=2, 4, 8, 16, 32, 64, or 128.
// They retain Phi for q > multiplier * tables.limit() and stream the bounded
// lower range. Alternatively, -DP2_HYBRID_SELECTOR (v1) or
// -DP2_HYBRID_SELECTOR_V2 or -DP2_HYBRID_SELECTOR_V3 chooses a frozen policy
// once per call. The four selection flags are mutually exclusive.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/ThreadBatch.h"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

#ifdef P2_HYBRID_MEMORY
#include "hybrid_metrics.h"
#endif

#if (defined(P2_HYBRID_SELECTOR) + defined(P2_HYBRID_SELECTOR_V2) + \
     defined(P2_HYBRID_SELECTOR_V3) + \
     defined(P2_HYBRID_MULTIPLIER)) > 1
#error "Select only one of P2_HYBRID_SELECTOR, P2_HYBRID_SELECTOR_V2, P2_HYBRID_SELECTOR_V3, or P2_HYBRID_MULTIPLIER"
#elif defined(P2_HYBRID_SELECTOR_V3)
#include "P2SelectorV3.h"
#elif defined(P2_HYBRID_SELECTOR_V2)
#include "P2SelectorV2.h"
#elif defined(P2_HYBRID_SELECTOR)
#include "P2Selector.h"
#elif !defined(P2_HYBRID_MULTIPLIER)
#error "Define P2_HYBRID_SELECTOR, P2_HYBRID_SELECTOR_V2, P2_HYBRID_SELECTOR_V3, or P2_HYBRID_MULTIPLIER as 2, 4, 8, 16, 32, 64, or 128"
#else
static_assert(P2_HYBRID_MULTIPLIER == 2 || P2_HYBRID_MULTIPLIER == 4 ||
              P2_HYBRID_MULTIPLIER == 8 || P2_HYBRID_MULTIPLIER == 16 ||
              P2_HYBRID_MULTIPLIER == 32 || P2_HYBRID_MULTIPLIER == 64 ||
              P2_HYBRID_MULTIPLIER == 128, "unsupported experimental multiplier");
#endif

namespace PrimeCounting {
namespace {
    using PartialSieveDetail::IntegerSqrt;

    std::int64_t countP2Range(std::int64_t x, std::int64_t first,
                              std::int64_t last, const PhiTables& tables) {
        const auto& primes = tables.primes();
        const auto& pi = tables.pi();
        std::int64_t sum = 0;
        for (std::int64_t i = first; i < last; ++i) {
            const std::int64_t quotient = x / primes[i];
            const std::int64_t c = pi[IntegerSqrt(quotient)];
            sum += Phi(quotient, c, tables) + c - 1 - (i - 1);
        }
        return sum;
    }

    std::int64_t countPhiPrefix(std::int64_t x, std::int64_t first,
                                std::int64_t last, std::size_t budget,
                                const PhiTables& tables) {
        if (first == last) return 0;
        const std::size_t workers = std::min<std::size_t>(
            budget, static_cast<std::size_t>(last - first));
#ifdef P2_HYBRID_MEMORY
        P2HybridMetrics::current().phiWorkers = workers;
#endif
        if (workers == 1) return countP2Range(x, first, last, tables);

        // Retain production's dynamic 32-query chunks. The budget derives
        // from the full original non-table range, so removing cheaper tail
        // queries does not by itself reduce workers on the costly prefix.
        constexpr std::int64_t CHUNK_SIZE = 32;
        std::atomic<std::int64_t> next(first);
        std::vector<std::int64_t> sums(workers, 0);
#ifdef P2_HYBRID_MEMORY
        P2HybridMetrics::current().phiSumBytes = sums.capacity() * sizeof(sums[0]);
#endif
        ThreadBatch batch(workers);
        for (std::size_t slot = 0; slot < workers; ++slot) {
            batch.launch([&, slot] {
                std::int64_t subtotal = 0;
                for (;;) {
                    const std::int64_t begin = next.fetch_add(CHUNK_SIZE,
                        std::memory_order_relaxed);
                    if (begin >= last) break;
                    subtotal += countP2Range(x, begin,
                        std::min(begin + CHUNK_SIZE, last), tables);
                }
                sums[slot] = subtotal;
            });
        }
        batch.finish();
        return std::accumulate(sums.begin(), sums.end(), std::int64_t(0));
    }

    // Adapted from the retained full-range streaming prototype. It counts
    // primes at increasing inclusive endpoints using one byte per odd integer
    // and persistent next-crossing positions for each required odd base prime.
    class PrimeCountStream {
        const std::vector<std::int64_t>& primes;
        std::vector<std::int64_t> next;
        std::vector<unsigned char> sieve;
        std::int64_t maximum, low, count;
        std::size_t offset = 0, length = 0;

        void fill() {
            low += 2 * static_cast<std::int64_t>(length);
            length = static_cast<std::size_t>(std::min<std::int64_t>(
                static_cast<std::int64_t>(sieve.size()), (maximum - low) / 2 + 1));
            std::fill(sieve.begin(), sieve.begin() + length, 1);
            const std::int64_t high = low + 2 * (static_cast<std::int64_t>(length) - 1);
            for (std::size_t i = 2; i < next.size(); ++i) {
                const std::int64_t p = primes[i];
                if (p > high / p) break;
                std::int64_t j = (next[i] - low) / 2;
                for (; j < static_cast<std::int64_t>(length); j += p) sieve[j] = 0;
                next[i] = low + 2 * j;
            }
            offset = 0;
        }

    public:
        PrimeCountStream(std::int64_t start, std::int64_t finish,
                         std::int64_t piStart, const PhiTables& tables)
            : primes(tables.primes()),
              next(static_cast<std::size_t>(tables.pi()[IntegerSqrt(finish)]) + 1),
              sieve(static_cast<std::size_t>(std::max<std::int64_t>(
                  32768, IntegerSqrt(finish)))),
              maximum(finish), low(start + 1), count(piStart) {
            if (low % 2 == 0) ++low;
            for (std::size_t i = 2; i < next.size(); ++i) {
                const std::int64_t p = primes[i];
                std::int64_t first = std::max(p * p, low + (p - low % p) % p);
                if (first % 2 == 0) first += p;
                next[i] = first;
            }
        }

        // All targets lie in [start, maximum] and are nondecreasing. A final
        // even target may leave low one past maximum without needing a fill.
        std::int64_t countTo(std::int64_t target) {
            while (low + 2 * static_cast<std::int64_t>(offset) <= target) {
                if (offset == length) fill();
                const std::size_t stop = static_cast<std::size_t>(
                    std::min<std::int64_t>(static_cast<std::int64_t>(length),
                                           (target - low) / 2 + 1));
                count = std::accumulate(sieve.begin() + offset,
                                        sieve.begin() + stop, count);
                offset = stop;
            }
            return count;
        }

#ifdef P2_HYBRID_MEMORY
        void storage(P2HybridMetrics::StreamStorage& metrics) const {
            metrics.sieveBytes = sieve.capacity() * sizeof(sieve[0]);
            metrics.nextBytes = next.capacity() * sizeof(next[0]);
        }
#endif
    };

    std::int64_t piAt(std::int64_t x, const PhiTables& tables) {
        if (x <= tables.limit()) return tables.pi()[x];
        const std::int64_t a = tables.pi()[IntegerSqrt(x)];
        return Phi(x, a, tables) + a - 1;
    }

    std::int64_t countStreamRange(std::int64_t x, std::int64_t first,
                                  std::int64_t last, std::int64_t start,
                                  const PhiTables& tables
#ifdef P2_HYBRID_MEMORY
                                  , P2HybridMetrics::StreamStorage& metrics
#endif
                                  ) {
        const auto& primes = tables.primes();
        const std::int64_t finish = x / primes[first];
        PrimeCountStream stream(start, finish, piAt(start, tables), tables);
#ifdef P2_HYBRID_MEMORY
        metrics.active = true;
        metrics.phiSeed = start > tables.limit();
        metrics.start = start;
        metrics.finish = finish;
        metrics.queries = last - first;
        stream.storage(metrics);
#endif
        std::int64_t sum = 0;
        for (std::int64_t i = last; i-- > first;) {
            sum += stream.countTo(x / primes[i]) - (i - 1);
        }
        return sum;
    }

    std::int64_t countStreamingTail(std::int64_t x, std::int64_t first,
                                    std::int64_t last, int nThreads,
                                    const PhiTables& tables) {
        if (first == last) return 0;
        const auto& primes = tables.primes();
        const std::int64_t span = x / primes[first] - tables.limit();
        constexpr std::int64_t MIN_WORK_PER_THREAD = 1048576;
        const int workers = static_cast<int>(std::min<std::int64_t>(
            std::max(1, nThreads), std::min(last - first,
                std::max<std::int64_t>(1, span / MIN_WORK_PER_THREAD))));
#ifdef P2_HYBRID_MEMORY
        auto& metrics = P2HybridMetrics::current();
        metrics.intervals.resize(workers);
#endif
        if (workers == 1) {
#ifdef P2_HYBRID_MEMORY
            metrics.streamWorkers = 1;
#endif
            return countStreamRange(x, first, last, tables.limit(), tables
#ifdef P2_HYBRID_MEMORY
                                    , metrics.intervals[0]
#endif
                                    );
        }

        // Each disjoint quotient interval is seeded once. Equal numeric spans
        // approximate the marking work better than equal prime-index counts.
        // The first interval uses pi(L); later seeds use serial Phi, never
        // nested workers. Empty intervals own no queries and launch no thread.
        const std::int64_t step = span / workers, extra = span % workers;
        std::int64_t lower = tables.limit();
        std::vector<std::int64_t> sums(workers, 0);
#ifdef P2_HYBRID_MEMORY
        metrics.streamSumBytes = sums.capacity() * sizeof(sums[0]);
#endif
        ThreadBatch batch(workers);
        for (int slot = 0; slot < workers && last > first; ++slot) {
            const std::int64_t upper = lower + step + (slot < extra);
            const std::int64_t begin = std::upper_bound(
                primes.begin() + first, primes.begin() + last, x / (upper + 1)
            ) - primes.begin();
            if (begin < last) {
#ifdef P2_HYBRID_MEMORY
                ++metrics.streamWorkers;
#endif
                batch.launch([&, slot, begin, last, lower] {
                    sums[slot] = countStreamRange(x, begin, last, lower, tables
#ifdef P2_HYBRID_MEMORY
                                                  , metrics.intervals[slot]
#endif
                                                  );
                });
            }
            last = begin;
            lower = upper;
        }
        batch.finish();
        return std::accumulate(sums.begin(), sums.end(), std::int64_t(0));
    }
} // namespace

    std::int64_t P2(std::int64_t x, std::int64_t a, int nThreads,
                    const PhiTables& tables) {
        PartialSieveDetail::Validate(x, a, tables);
#ifdef P2_HYBRID_MEMORY
        auto& metrics = P2HybridMetrics::current();
        metrics = P2HybridMetrics::Snapshot();
        metrics.x = x;
        metrics.a = a;
        metrics.tableLimit = tables.limit();
#endif
        if (x < 4) return 0;
        const std::int64_t root = IntegerSqrt(x);
        if (tables.limit() < root) {
            throw std::out_of_range("P2 requires tables through floor(sqrt(x))");
        }
        const auto& primes = tables.primes();
        const auto& pi = tables.pi();
        const std::int64_t b = pi[root];
        if (a >= b) return 0;
        const std::int64_t first = a + 1;
        const std::int64_t last = std::upper_bound(
            primes.begin() + first, primes.begin() + b + 1,
            x / (tables.limit() + 1)) - primes.begin();

        std::int64_t sum = 0;
        for (std::int64_t i = last; i <= b; ++i) {
            sum += static_cast<std::int64_t>(pi[x / primes[i]]) - (i - 1);
        }
#ifdef P2_HYBRID_MEMORY
        metrics.directQueries = b + 1 - last;
#endif
        if (last == first) return sum;

#if defined(P2_HYBRID_SELECTOR) || defined(P2_HYBRID_SELECTOR_V2) || \
    defined(P2_HYBRID_SELECTOR_V3)
        constexpr std::int64_t MIN_TERMS_PER_THREAD = 2048;
        const std::size_t phiBudget = static_cast<std::size_t>(
            std::min<std::int64_t>(std::max(1, nThreads),
                std::max<std::int64_t>(1, (last - first) / MIN_TERMS_PER_THREAD)));
#ifdef P2_HYBRID_SELECTOR_V3
        const P2SelectorV3::Decision choice = P2SelectorV3::choose(
            x, tables.limit(), first, last, phiBudget, nThreads, primes);
#elif defined(P2_HYBRID_SELECTOR_V2)
        const P2SelectorV2::Decision choice = P2SelectorV2::choose(
            x, tables.limit(), first, last, phiBudget, nThreads, primes);
#else
        const P2Selector::Decision choice = P2Selector::choose(
            x, tables.limit(), first, last, phiBudget, primes);
#endif
        const std::int64_t streamFirst = choice.streamFirst;
#else
        const std::int64_t largestQuotient = x / primes[first];
        // Cap before multiplying, so forming Q cannot overflow even if the
        // accepted table range is expanded in a future experiment.
        const std::int64_t streamLimit =
            tables.limit() > largestQuotient / P2_HYBRID_MULTIPLIER ?
            largestQuotient : tables.limit() * P2_HYBRID_MULTIPLIER;
        const std::int64_t streamFirst = std::upper_bound(
            primes.begin() + first, primes.begin() + last,
            x / (streamLimit + 1)) - primes.begin();
#endif
#ifdef P2_HYBRID_MEMORY
#if defined(P2_HYBRID_SELECTOR) || defined(P2_HYBRID_SELECTOR_V2) || \
    defined(P2_HYBRID_SELECTOR_V3)
        metrics.streamLimit = choice.limit;
#else
        metrics.streamLimit = streamLimit;
#endif
        metrics.phiQueries = streamFirst - first;
        metrics.streamQueries = last - streamFirst;
#endif

#if !defined(P2_HYBRID_SELECTOR) && !defined(P2_HYBRID_SELECTOR_V2) && \
    !defined(P2_HYBRID_SELECTOR_V3)
        constexpr std::int64_t MIN_TERMS_PER_THREAD = 2048;
        const std::size_t phiBudget = static_cast<std::size_t>(
            std::min<std::int64_t>(std::max(1, nThreads),
                std::max<std::int64_t>(1, (last - first) / MIN_TERMS_PER_THREAD)));
#endif
        sum += countPhiPrefix(x, first, streamFirst, phiBudget, tables);
        // Prefix workers and storage are finished before tail workers start.
        // Both phases borrow the same immutable context and process-wide cache.
        sum += countStreamingTail(x, streamFirst, last, nThreads, tables);
        return sum;
    }
}
