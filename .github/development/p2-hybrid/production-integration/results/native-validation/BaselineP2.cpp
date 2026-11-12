#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/ThreadBatch.h"
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <numeric>
#include <vector>

namespace PrimeCounting {
namespace {

    std::int64_t countP2Range(std::int64_t x, std::int64_t first,
                              std::int64_t last, const PhiTables& tables) {
        const auto& primes = tables.primes();
        const auto& pi = tables.pi();
        std::int64_t sum = 0;

        for (std::int64_t i = first; i < last; ++i) {
            const std::int64_t quotient = x / primes[i];
            const std::int64_t c = pi[PartialSieveDetail::IntegerSqrt(quotient)];
            // q >= p[i] >= 2. Legendre's identity gives pi(q) using the
            // existing serial Phi evaluator and the caller's larger table.
            // isqrt(q) <= isqrt(x), so its cutoff is always prepared.
            const std::int64_t piQuotient = Phi(quotient, c, tables) + c - 1;
            sum += piQuotient - (i - 1);
        }

        return sum;
    }
}

    std::int64_t P2(std::int64_t x, std::int64_t a, int nThreads,
                    const PhiTables& tables) {
        PartialSieveDetail::Validate(x, a, tables);
        if (x < 4) return 0;

        const std::int64_t root = PartialSieveDetail::IntegerSqrt(x);

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
            x / (tables.limit() + 1)
        ) - primes.begin();

        std::int64_t sum = 0;

        for (std::int64_t i = last; i <= b; ++i) {
            sum += static_cast<std::int64_t>(pi[x / primes[i]]) - (i - 1);
        }

        if (last == first) return sum;

        // Tiny calls cost less than launching workers. This is an empirical
        // grain, not a restriction on the mathematical input domain.
        constexpr std::int64_t MIN_TERMS_PER_THREAD = 2048;

        const std::size_t workers = static_cast<std::size_t>(
            std::min<std::int64_t>(std::max(1, nThreads),
                std::max<std::int64_t>(1, (last - first) / MIN_TERMS_PER_THREAD))
        );

        if (workers == 1) return sum + countP2Range(x, first, last, tables);

        // Early quotients can be much more expensive than later ones. Workers
        // take small contiguous chunks until the sum is exhausted. The atomic
        // counter assigns unique indices only; it publishes no other data, so
        // relaxed ordering suffices. Launch/join synchronize tables/results.
        constexpr std::int64_t CHUNK_SIZE = 32;
        std::atomic<std::int64_t> next(first);
        std::vector<std::int64_t> sums(workers, 0);
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
        return sum + std::accumulate(sums.begin(), sums.end(), std::int64_t(0));
    }
}
