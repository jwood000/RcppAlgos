#include "NumbersUtils/Phi.h"
#include "NumbersUtils/Eratosthenes.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <mutex>
#include <new>
#include <stdexcept>

// General partial-sieve evaluator developed from PiLegendre.cpp's
// recurrence and shortcuts (based on Kim Walisch's primecount implementation).
// PhiParallel partitions the outer sum; each worker uses this serial evaluator.

namespace PrimeCounting {
namespace {

    constexpr std::int64_t MAX_PHI_X = (std::int64_t(1) << 53) - 1;
    // trunc(sqrt(2^53 - 1)) = 94906265
    constexpr std::int64_t MAX_TABLE_LIMIT = 94906265;
    constexpr int TINY_A = 6;
    constexpr std::array<int, 7> primeProducts = {
        {1, 2, 6, 30, 210, 2310, 30030}
    };
    constexpr std::array<int, 7> totients = {
        {1, 1, 2, 8, 48, 480, 5760}
    };

    std::int64_t integerSqrt(std::int64_t x) {
        if (x == 0) return 0;
        std::int64_t root = static_cast<std::int64_t>(
            std::sqrt(static_cast<double>(x))
        );
        // Correct floating-point rounding at square boundaries using integers.
        while (root > x / root) --root;
        while (root + 1 <= x / (root + 1)) ++root;
        return root;
    }

    // Phi(x, a) counts integers in [1, x] not divisible by the first a
    // primes. For a >= 1, divisibility repeats every P = primeProducts[a],
    // and each full period contributes totients[a] survivors.
    // counts[a][r] stores the number of survivors in the prefix [1, r],
    // for 0 <= r < P, so a partial period needs only one lookup.
    struct TinyPhiTables {
        std::array<std::vector<std::uint16_t>, TINY_A + 1> counts;

        TinyPhiTables() {
            constexpr std::array<int, TINY_A> primes = {{2, 3, 5, 7, 11, 13}};
            for (int a = 0; a <= TINY_A; ++a) {
                const int period = primeProducts[a];

                // This is a partial sieve, not a primality test: 1 survives,
                // as do composites with no factor among these a primes
                // (e.g. 121 survives when a = 4 removes 2, 3, 5, and 7).
                std::vector<char> survives(period, 1);
                survives[0] = 0;

                for (int i = 0; i < a; ++i) {
                    for (int n = primes[i]; n < period; n += primes[i]) {
                        survives[n] = 0;
                    }
                }

                // Convert individual survivor flags into cumulative counts.
                // counts[a][0] = 0 represents the empty prefix. For a = 2,
                // survives is {0,1,0,0,0,1}, but counts is {0,1,1,1,1,2}.
                counts[a].resize(period, 0);

                for (int n = 1; n < period; ++n) {
                    counts[a][n] = counts[a][n - 1] + survives[n];
                }
            }
        }
    };

    std::int64_t phiTiny(std::int64_t x, std::int64_t a) {
        if (a == 0) return x;

        // One shared object, initialized on the first call reaching here,
        // not once per call or thread. C++11 synchronizes initialization;
        // all later evaluations only read the immutable tables.
        static const TinyPhiTables tables;
        const std::int64_t period = primeProducts[a];

        // Count complete periods, then the survivors in the remaining prefix.
        // For x = 14, a = 2: 2 * 2 + counts[2][2] = 5.
        return (x / period) * totients[a] + tables.counts[a][x % period];
    }

    // Bounded cache of positive Phi(x, a) answers, shared across evaluations:
    //   PhiCache -> 100 rows indexed by a
    //   PhiCacheRow -> 16 atomic page pointers, indexed by x / 4096
    //   PhiPage -> 4096 contiguous atomic cells, indexed by x % 4096
    // Pages are allocated only when needed; out-of-range inputs are still
    // evaluated normally, without caching. Zero denotes a cache miss.
    // For positive x, 1 <= Phi(x, a) <= x, so cached answers fit in uint16_t.
    // Public Phi handles x = 0 before consulting or populating the cache.
    constexpr std::size_t CACHE_ROWS = 100;
    constexpr unsigned PAGE_BITS = 12;
    constexpr std::size_t PAGE_SIZE = 1u << PAGE_BITS;
    constexpr std::size_t PAGE_COUNT = 65536u / PAGE_SIZE;

    struct PhiCell {
        std::atomic<std::uint16_t> value;
        // The wrapper guarantees an initial cache miss for every array cell.
        // A default-constructed C++11 atomic does not initialize its value.
        PhiCell() : value(0) {}
    };

    struct PhiPage {
        // One allocation holds all cells inline; no allocation per PhiCell.
        PhiCell cells[PAGE_SIZE];
    };

    struct PublishedPage {
        // Atomic access protects the shared pointer, not the pointed-to data.
        // PhiCell provides the separate atomic access needed for each answer.
        std::atomic<PhiPage*> value;
        PublishedPage() : value(nullptr) {}
    };

    struct PhiCacheRow {
        // Serializes page creation within this row only. Readers never take
        // this mutex, and writers skip it when their page already exists.
        std::mutex allocationMutex;
        PublishedPage pages[PAGE_COUNT];

        ~PhiCacheRow() {
            // Published pages never move or get freed during evaluation.
            // Destruction requires all evaluations to have finished; an
            // atomic pointer alone would not protect readers from deletion.
            // With no concurrent users, relaxed loads suffice for cleanup.
            for (auto& page : pages) {
                delete page.value.load(std::memory_order_relaxed);
            }
        }
    };

    class PhiCache {
        std::array<PhiCacheRow, CACHE_ROWS> rows;

    public:
        std::int64_t get(std::uint64_t x, std::uint64_t a) const {
            if (a >= CACHE_ROWS || x > std::numeric_limits<std::uint16_t>::max()) {
                return 0;
            }

            // The shift and mask implement division and remainder by 4096.
            // Acquiring a published pointer makes its page initialization
            // visible, even though this reader never takes the row mutex.
            const PhiPage* page = rows[a].pages[x >> PAGE_BITS].value.load(
                std::memory_order_acquire
            );

            // The cell is the complete answer, not a ready flag for other
            // data. Relaxed access is still atomic. Seeing zero simply makes
            // the caller compute the answer, even if another worker has it.
            return page ? page->cells[x & (PAGE_SIZE - 1)].value.load(
                std::memory_order_relaxed
            ) : 0;
        }

        void put(std::uint64_t x, std::uint64_t a, std::int64_t value) {
            if (a >= CACHE_ROWS || x > std::numeric_limits<std::uint16_t>::max()) {
                return;
            }

            auto& row = rows[a];
            auto& published = row.pages[x >> PAGE_BITS].value;
            PhiPage* page = published.load(std::memory_order_acquire);

            if (!page) {
                // lock_guard releases the mutex on every exit from this scope,
                // including an early return or a C++ exception.
                std::lock_guard<std::mutex> lock(row.allocationMutex);

                // Another writer may have created the page while we waited.
                page = published.load(std::memory_order_acquire);

                if (!page) {
                    page = new (std::nothrow) PhiPage;
                    if (!page) return; // Caching is optional.

                    // Construct and zero all cells before sharing the address.
                    // A reader's acquire load of this released pointer sees
                    // that initialization. The mutex alone cannot do this for
                    // readers, since they do not lock it.
                    published.store(page, std::memory_order_release);
                }
            }

            // A cell holds the complete positive answer. Duplicate writers
            // agree, even when their prepared tables have different limits.
            // Concurrent writes still require atomics, even for equal values.
            // No other data is published by this store, so relaxed is enough.
            // Duplicate computation is allowed; workers need not wait for an
            // in-progress answer. The sign belongs to the caller's recurrence.
            page->cells[x & (PAGE_SIZE - 1)].value.store(
                static_cast<std::uint16_t>(std::abs(value)),
                std::memory_order_relaxed
            );
        }
    };

    PhiCache& sharedCache() {
        // One persistent cache for all callers and workers. C++11 synchronizes
        // its first initialization; its mutable contents use the rules above.
        static PhiCache cache;
        return cache;
    }

    class PhiEvaluator {
        const std::vector<std::int64_t>& primes;
        const std::vector<std::int32_t>& pi;
        PhiCache& cache;

    public:
        PhiEvaluator(const std::vector<std::int64_t>& preparedPrimes,
                     const std::vector<std::int32_t>& preparedPi)
            : primes(preparedPrimes), pi(preparedPi), cache(sharedCache()) {}

        // Private kernel: a >= 6, x > p[a], and the direct pi shortcut has
        // already been ruled out. Public Phi checks the initial call.
        // A recursive call exists only when b > min(b, 6), so i >= 6.
        // Also p[i + 1] <= floor(sqrt(x)), hence its quotient is at least
        // p[i + 1] > p[i]. The parent tests the shortcut before recursing.
        // These invariants let children skip the initial checks safely.
        template <std::int64_t SIGN>
        std::int64_t evaluate(std::int64_t x, std::int64_t a) {
            if (a == TINY_A) return phiTiny(x, TINY_A) * SIGN;

            const std::int64_t cached = cache.get(x, a);
            if (cached) return cached * SIGN;

            const std::int64_t root = integerSqrt(x);
            std::int64_t b = a;

            if (root < static_cast<std::int64_t>(pi.size())) {
                b = std::min(a, static_cast<std::int64_t>(pi[root]));
            }

            // Otherwise every prepared prime is below root, so b = a exactly.
            // Keep the general base here: b can be below 6 even though a >= 6.
            // For example, Phi(18, 7) with a table through 17 has b = c = 2.
            const std::int64_t c = std::min(
                b, static_cast<std::int64_t>(TINY_A)
            );

            // x > p[a] here. Each term with b <= i < a equals one because
            // sqrt(x) < p[i + 1] <= p[a] < x. No zero-quotient tail is omitted.
            std::int64_t sum = (b - a + phiTiny(x, c)) * SIGN;
            std::int64_t i = c;

            for (; i < b; ++i) {
                const std::int64_t quotient = x / primes[i + 1];
                // i + 1 <= b <= a guarantees the lookahead prime exists.
                // A child is called only after this shortcut has failed.
                if (quotient < static_cast<std::int64_t>(pi.size()) &&
                    quotient < primes[i + 1] * primes[i + 1]) break;
                sum += evaluate<-SIGN>(quotient, i);
            }

            // Every i < b has a prepared p[i + 1]. Quotients decrease and
            // prime squares increase, so the direct shortcut stays valid.
            for (; i < b; ++i) {
                const std::int64_t quotient = x / primes[i + 1];
                sum += (static_cast<std::int64_t>(pi[quotient]) - i + 1) * -SIGN;
            }

            cache.put(x, a, sum);
            return sum;
        }
    };
}

    PhiTables::PhiTables(std::int64_t limit) {
        if (limit < 0 || limit > MAX_TABLE_LIMIT) {
            throw std::out_of_range("PhiTables limit must be between 0 and 94906265");
        }

        PrimeSieve::sqrtBigPrimes(static_cast<int>(limit), true, false, true, primes_);
        pi_.resize(static_cast<std::size_t>(limit) + 1, 0);
        std::int32_t count = 0;

        for (std::int64_t n = 1; n <= limit; ++n) {
            if (static_cast<std::size_t>(count + 1) < primes_.size() &&
                primes_[count + 1] == n) {
                ++count;
            }
            pi_[n] = count;
        }
    }

    std::int64_t PhiTables::limit() const {
        return static_cast<std::int64_t>(pi_.size()) - 1;
    }

    std::int64_t PhiTables::primeCount() const {
        return static_cast<std::int64_t>(primes_.size()) - 1;
    }

    std::int64_t Phi(std::int64_t x, std::int64_t a, const PhiTables& tables) {
        if (x < 0 || x > MAX_PHI_X) {
            throw std::out_of_range("Phi x must be between 0 and 2^53 - 1");
        }

        if (a < 0 || a > tables.primeCount()) {
            throw std::out_of_range("Phi a must be between 0 and the prepared prime count");
        }

        if (x == 0) return 0;
        if (a == 0) return x;
        if (x <= tables.primes_[a]) return 1;
        if (a <= TINY_A) return phiTiny(x, a);

        // Check the initial direct shortcut here; recursive children are
        // checked by their parent. A final prepared prime may have no
        // lookahead prime, which disables only this shortcut, not evaluation.
        if (x < static_cast<std::int64_t>(tables.pi_.size()) &&
            a + 1 < static_cast<std::int64_t>(tables.primes_.size()) &&
            x < tables.primes_[a + 1] * tables.primes_[a + 1]) {
            return static_cast<std::int64_t>(tables.pi_[x]) - a + 1;
        }

        PhiEvaluator evaluator(tables.primes_, tables.pi_);
        return evaluator.evaluate<1>(x, a);
    }
}
