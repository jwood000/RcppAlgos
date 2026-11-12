#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/ThreadBatch.h"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace PrimeCounting {
namespace {
    using PartialSieveDetail::IntegerSqrt;

    // Selection and launching must use the same numeric-span worker grain.
    constexpr std::int64_t MIN_STREAM_WORK_PER_THREAD = 1048576;

    // The earliest Phi queries can be much more expensive than the rest.
    // Give the first 32 queries individual tickets, then use 32-query jobs.
    // One queue covers both parts: workers can claim later jobs while an
    // earlier query is still running. No extra barrier or worker is needed.
    class PhiWorkQueue {
        const std::int64_t first_, last_, headTerms_;
        const std::uint64_t headTasks_, tasks_;
        std::atomic<std::uint64_t> next_;

        static std::int64_t length(std::int64_t first, std::int64_t last) {
            if (first < 0 || last < first)
                throw std::out_of_range("invalid P2 work interval");
            return last - first;
        }

        static std::uint64_t tasks(std::int64_t terms, std::int64_t grain) {
            // Avoid overflowing terms+grain-1 when rounding up.
            return static_cast<std::uint64_t>(terms / grain + (terms % grain != 0));
        }

    public:
        PhiWorkQueue(std::int64_t first, std::int64_t last)
            : first_(first), last_(last),
              headTerms_(std::min<std::int64_t>(32, length(first, last))),
              headTasks_(static_cast<std::uint64_t>(headTerms_)),
              tasks_(headTasks_ + tasks(last - first - headTerms_, 32)), next_(0) {}

        std::uint64_t size() const { return tasks_; }

        // Disjoint ticket ranges cover [first,last) exactly once. Check the
        // ticket before multiplying; every valid offset is below last-first,
        // so adding first is safe even for a range ending at INT64_MAX.
        bool range(std::uint64_t ticket, std::int64_t& begin,
                   std::int64_t& end) const {
            if (ticket >= tasks_) return false;
            if (ticket < headTasks_) {
                const auto offset = static_cast<std::int64_t>(ticket);
                begin = first_ + offset;
                end = begin + 1;
            } else {
                const auto offset = headTerms_ +
                    static_cast<std::int64_t>(ticket - headTasks_) * 32;
                begin = first_ + offset;
                end = begin + std::min<std::int64_t>(32, last_ - begin);
            }
            return true;
        }

        bool claim(std::int64_t& begin, std::int64_t& end) {
            // The counter assigns work only. ThreadBatch launch/join handles
            // publication and result visibility, so relaxed ordering suffices.
            return range(next_.fetch_add(1, std::memory_order_relaxed), begin, end);
        }
    };

    // Choose a bounded tail of non-table quotients for segmented counting.
    // Returning last keeps the entire range on Phi. These empirical guards
    // use only the requested work, never method labels, cache state or timings.
    // Input/table validation has already established one-based prime indices.
    std::int64_t chooseStreamFirst(std::int64_t x, std::int64_t limit,
                                    std::int64_t first, std::int64_t last,
                                    std::size_t phiBudget, int nThreads,
                                    const std::vector<std::int64_t>& primes) {
        if (last - first <= 1) return last;

        int multiplier;
        if (x < 100000000000LL) {
            multiplier = 8;
        } else if (x < 3000000000000LL) {
            multiplier = 16;
        } else if (x < 30000000000000LL) {
            multiplier = 32;
        } else {
            multiplier = 64;
        }

        if (multiplier == 16 && x >= 1000000000000LL && phiBudget >= 4) {
            multiplier = 32;
        }

        const bool smallParallel = x < 100000000000LL && phiBudget > 1;
        if (smallParallel) multiplier = 16;

        const std::int64_t largestQuotient = x / primes[first];
        // Cap before multiplying. floor(x/p) <= Q is equivalent to
        // p > floor(x/(Q+1)); Q <= largestQuotient <= x/2, so Q+1 is safe.
        const std::int64_t streamLimit = limit > largestQuotient / multiplier ?
            largestQuotient : limit * multiplier;
        const std::int64_t streamFirst = std::upper_bound(
            primes.begin() + first, primes.begin() + last,
            x / (streamLimit + 1)) - primes.begin();
        const std::int64_t streamQueries = last - streamFirst;
        std::int64_t span = 0;
        std::size_t streamBudget = 0;

        if (streamQueries > 0) {
            // Use the actual final target, which can be below the cap. Match
            // the launcher's span grain; empty bands may use fewer workers.
            span = x / primes[streamFirst] - limit;
            streamBudget = static_cast<std::size_t>(
                std::min<std::int64_t>(std::max(1, nThreads),
                    std::min(streamQueries,
                        std::max<std::int64_t>(1, span / MIN_STREAM_WORK_PER_THREAD))));
        }

        // A short complete range needs no Phi prefix and at most two standard
        // 32768-odd-number segments. This exception precedes the other guards.
        if (streamFirst == first && span <= 131072) return streamFirst;
        if (last - first < 4096) return last;
        if (streamQueries < 4096) return last;
        if (smallParallel && streamBudget < phiBudget) return last;

        if (smallParallel) {
            const auto root = IntegerSqrt(x);
            if (limit != root) return last;
            // Conservative admission: root-sized tables with a cutoff at
            // least as large as Lehmer's. Earlier cutoffs and larger tables
            // have less consistent first/warmed benefits from streaming.
            // This is an empirical limit, not a mathematical crossover.
            if (primes[first] <= IntegerSqrt(root)) return last;
        }
        return streamFirst;
    }

    std::int64_t countP2Range(std::int64_t x, std::int64_t first,
                              std::int64_t last, const PhiTables& tables) {
        const auto& primes = tables.primes();
        const auto& pi = tables.pi();
        std::int64_t sum = 0;
        for (std::int64_t i = first; i < last; ++i) {
            const std::int64_t quotient = x / primes[i];
            const std::int64_t c = pi[IntegerSqrt(quotient)];
            // q >= p[i] >= 2. Legendre's identity gives pi(q) using the
            // serial Phi evaluator and the caller's larger table.
            // isqrt(q) <= isqrt(x), so its cutoff is always prepared.
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
        if (workers == 1) return countP2Range(x, first, last, tables);

        // Phi-only calls and hybrid prefixes share the same initial jobs.
        // The budget derives from the full original non-table range, so
        // removing cheaper tail queries does not reduce the prefix's budget.
        PhiWorkQueue next(first, last);
        std::vector<std::int64_t> sums(workers, 0);
        ThreadBatch batch(workers);
        for (std::size_t slot = 0; slot < workers; ++slot) {
            batch.launch([&, slot] {
                std::int64_t subtotal = 0;
                for (;;) {
                    std::int64_t begin, end;
                    if (!next.claim(begin, end)) break;
                    subtotal += countP2Range(x, begin, end, tables);
                }
                sums[slot] = subtotal;
            });
        }
        batch.finish();
        return std::accumulate(sums.begin(), sums.end(), std::int64_t(0));
    }

    // Count primes at increasing inclusive endpoints using one byte per odd
    // integer and persistent next-crossing positions for each odd base prime.
    // Each worker owns its stream; only the prepared prime vector is shared.
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
    };

    std::int64_t piAt(std::int64_t x, const PhiTables& tables) {
        if (x <= tables.limit()) return tables.pi()[x];
        const std::int64_t a = tables.pi()[IntegerSqrt(x)];
        return Phi(x, a, tables) + a - 1;
    }

    std::int64_t countStreamRange(std::int64_t x, std::int64_t first,
                                  std::int64_t last, std::int64_t start,
                                  const PhiTables& tables) {
        const auto& primes = tables.primes();
        const std::int64_t finish = x / primes[first];
        PrimeCountStream stream(start, finish, piAt(start, tables), tables);
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
        const int workers = static_cast<int>(std::min<std::int64_t>(
            std::max(1, nThreads), std::min(last - first,
                std::max<std::int64_t>(1, span / MIN_STREAM_WORK_PER_THREAD))));
        if (workers == 1) {
            return countStreamRange(x, first, last, tables.limit(), tables);
        }

        // Each disjoint quotient interval is seeded once. Equal numeric spans
        // approximate the marking work better than equal prime-index counts.
        // The first interval uses pi(L); later seeds use serial Phi, never
        // nested workers. Empty intervals own no queries and launch no thread.
        const std::int64_t step = span / workers, extra = span % workers;
        std::int64_t lower = tables.limit();
        std::vector<std::int64_t> sums(workers, 0);
        ThreadBatch batch(workers);
        for (int slot = 0; slot < workers && last > first; ++slot) {
            const std::int64_t upper = lower + step + (slot < extra);
            const std::int64_t begin = std::upper_bound(
                primes.begin() + first, primes.begin() + last, x / (upper + 1)
            ) - primes.begin();
            if (begin < last) {
                batch.launch([&, slot, begin, last, lower] {
                    sums[slot] = countStreamRange(x, begin, last, lower, tables);
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

        // Each nondecreasing pair p[i] <= p[j] counts once, including squares:
        // P2 = sum(i=a+1..b, pi(floor(x/p[i])) - (i-1)). Quotients decrease
        // with i; separate the direct-table suffix once, outside the hot loop.
        // floor(x/p) <= limit exactly when p > floor(x/(limit+1)).
        const std::int64_t last = std::upper_bound(
            primes.begin() + first, primes.begin() + b + 1,
            x / (tables.limit() + 1)) - primes.begin();

        std::int64_t sum = 0;
        for (std::int64_t i = last; i <= b; ++i) {
            sum += static_cast<std::int64_t>(pi[x / primes[i]]) - (i - 1);
        }
        if (last == first) return sum;

        // Base the Phi budget on the full non-table range. Removing cheaper
        // tail queries must not by itself reduce workers on the costly prefix.
        constexpr std::int64_t MIN_TERMS_PER_THREAD = 2048;
        const std::size_t phiBudget = static_cast<std::size_t>(
            std::min<std::int64_t>(std::max(1, nThreads),
                std::max<std::int64_t>(1, (last - first) / MIN_TERMS_PER_THREAD)));
        const std::int64_t streamFirst = chooseStreamFirst(
            x, tables.limit(), first, last, phiBudget, nThreads, primes);

        sum += countPhiPrefix(x, first, streamFirst, phiBudget, tables);
        // Prefix workers and storage finish before tail workers start. Both
        // phases borrow the immutable context and process-wide atomic cache.
        sum += countStreamingTail(x, streamFirst, last, nThreads, tables);
        return sum;
    }
}
