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
// Legendre uses this evaluator while retaining its own work scheduling.

namespace PrimeCounting {
namespace {

    constexpr std::int64_t MAX_PHI_X = (std::int64_t(1) << 53) - 1;
    constexpr std::int64_t MAX_TABLE_LIMIT = 94906265;
    constexpr int TINY_A = 6;
    constexpr std::array<int, 7> primeProducts = {{1, 2, 6, 30, 210, 2310, 30030}};
    constexpr std::array<int, 7> totients = {{1, 1, 2, 8, 48, 480, 5760}};

    std::int64_t integerSqrt(std::int64_t x) {
        if (x == 0) return 0;
        std::int64_t root = static_cast<std::int64_t>(std::sqrt(static_cast<double>(x)));
        // Correct floating-point rounding at square boundaries using integers.
        while (root > x / root) --root;
        while (root + 1 <= x / (root + 1)) ++root;
        return root;
    }

    struct TinyPhiTables {
        std::array<std::vector<std::uint16_t>, TINY_A + 1> counts;

        TinyPhiTables() {
            constexpr std::array<int, TINY_A> primes = {{2, 3, 5, 7, 11, 13}};
            for (int a = 0; a <= TINY_A; ++a) {
                const int period = primeProducts[a];
                std::vector<char> survives(period, 1);
                survives[0] = 0;
                for (int i = 0; i < a; ++i) {
                    for (int n = primes[i]; n < period; n += primes[i]) {
                        survives[n] = 0;
                    }
                }
                counts[a].resize(period, 0);
                for (int n = 1; n < period; ++n) {
                    counts[a][n] = counts[a][n - 1] + survives[n];
                }
            }
        }
    };

    std::int64_t phiTiny(std::int64_t x, std::int64_t a) {
        if (a == 0) return x;
        // Generate the small periodic tables lazily. Initialization is
        // thread-safe under C++11; subsequent evaluations only read them.
        static const TinyPhiTables tables;
        const std::int64_t period = primeProducts[a];
        return (x / period) * totients[a] + tables.counts[a][x % period];
    }

    constexpr std::size_t CACHE_ROWS = 100;
    constexpr unsigned PAGE_BITS = 12;
    constexpr std::size_t PAGE_SIZE = 1u << PAGE_BITS;
    constexpr std::size_t PAGE_COUNT = 65536u / PAGE_SIZE;

    struct PhiCell {
        std::atomic<std::uint16_t> value;
        PhiCell() : value(0) {}
    };

    struct PhiPage {
        PhiCell cells[PAGE_SIZE];
    };

    struct PublishedPage {
        std::atomic<PhiPage*> value;
        PublishedPage() : value(nullptr) {}
    };

    struct PhiCacheRow {
        std::mutex allocationMutex;
        PublishedPage pages[PAGE_COUNT];

        ~PhiCacheRow() {
            // Destruction requires all evaluations to have finished.
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
            const PhiPage* page = rows[a].pages[x >> PAGE_BITS].value.load(
                std::memory_order_acquire
            );
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
                std::lock_guard<std::mutex> lock(row.allocationMutex);
                page = published.load(std::memory_order_acquire);
                if (!page) {
                    page = new (std::nothrow) PhiPage;
                    if (!page) return; // Caching is optional.
                    published.store(page, std::memory_order_release);
                }
            }
            // A cell holds the complete positive answer. Duplicate writers
            // agree, even when their prepared tables have different limits.
            page->cells[x & (PAGE_SIZE - 1)].value.store(
                static_cast<std::uint16_t>(std::abs(value)),
                std::memory_order_relaxed
            );
        }
    };

    PhiCache& sharedCache() {
        static PhiCache cache;
        return cache;
    }

    class PhiEvaluator {
        const std::vector<std::int64_t>& primes;
        const std::vector<std::int32_t>& pi;
        PhiCache& cache;

        bool isPi(std::int64_t x, std::int64_t a) const {
            // A final prepared prime need not have a lookahead prime. Missing
            // lookahead disables this shortcut, not the recursive evaluation.
            return x < static_cast<std::int64_t>(pi.size()) &&
                   a + 1 < static_cast<std::int64_t>(primes.size()) &&
                   x < primes[a + 1] * primes[a + 1];
        }

    public:
        PhiEvaluator(const std::vector<std::int64_t>& preparedPrimes,
                     const std::vector<std::int32_t>& preparedPi)
            : primes(preparedPrimes), pi(preparedPi), cache(sharedCache()) {}

        template <std::int64_t SIGN>
        std::int64_t evaluate(std::int64_t x, std::int64_t a) {
            if (x == 0) return 0;
            if (a == 0) return x * SIGN;
            if (x <= primes[a]) return SIGN;
            if (a <= TINY_A) return phiTiny(x, a) * SIGN;
            if (isPi(x, a)) {
                return (static_cast<std::int64_t>(pi[x]) - a + 1) * SIGN;
            }

            const std::int64_t cached = cache.get(x, a);
            if (cached) return cached * SIGN;

            const std::int64_t root = integerSqrt(x);
            std::int64_t b = a;
            if (root < static_cast<std::int64_t>(pi.size())) {
                b = std::min(a, static_cast<std::int64_t>(pi[root]));
            }
            // Otherwise every prepared prime is below root, so b = a exactly.
            const std::int64_t c = std::min(b, static_cast<std::int64_t>(TINY_A));

            // x > p[a] here. Each term with b <= i < a equals one because
            // sqrt(x) < p[i + 1] <= p[a] < x. No zero-quotient tail is omitted.
            std::int64_t sum = (b - a + phiTiny(x, c)) * SIGN;
            std::int64_t i = c;
            for (; i < b; ++i) {
                const std::int64_t quotient = x / primes[i + 1];
                if (isPi(quotient, i)) break;
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

        PhiEvaluator evaluator(tables.primes_, tables.pi_);
        return evaluator.evaluate<1>(x, a);
    }
}
