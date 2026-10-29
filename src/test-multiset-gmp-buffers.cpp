#include "Partitions/BigPartsCountMultiset.h"
#include "Partitions/BigPartsCountSection.h"
#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>
#include <testthat.h>

static unsigned long Choose(int n, int k) {
    unsigned long result = 1;
    for (int i = 1; i <= k; ++i) result = result * (n - i + 1) / i;
    return result;
}

// Enumerate multiplicities independently of either DP recurrence.
static unsigned long BruteCount(
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& reps, std::size_t i, bool composition
) {
    if (i == allowed.size()) return n == 0 && m == 0;
    const int limit = std::max(0, std::min({reps[i], m, n / allowed[i]}));
    unsigned long result = 0;
    for (int t = 0; t <= limit; ++t) {
        result += BruteCount(n - t * allowed[i], m - t, allowed,
                             reps, i + 1, composition) *
            (composition ? Choose(m, t) : 1);
    }
    return result;
}

context("GMP multiset buffers") {
    test_that("dirty buffers are reset and storage is reused") {
        using Counter = void (*)(
            mpz_class&, std::vector<mpz_class>&, std::vector<mpz_class>&,
            int, int, const std::vector<int>&, const std::vector<int>&, int
        );
        const Counter counters[] = {
            CountPartsMultiset, CountCompsMultiset,
            CountCompsMultisetZero, CountCompsMultisetWeak
        };
        const std::size_t size = 41 * 9;
        mpz_class poison = 1;
        poison <<= 128;
        std::vector<mpz_class> p1(size, poison), p2(size, -poison);
        const auto* first = p1.data();
        const auto* second = p2.data();
        std::mt19937 rng(104729);

        for (int trial = 0; trial < 5000; ++trial) {
            const int n = rng() % 41;
            const int m = rng() % 9;
            const int k = rng() % 9;
            const int strtLen = rng() % (m + 2);
            std::vector<int> allowed(k), reps(k + 2);
            int value = 0;
            for (int& a : allowed) {
                value += 1 + rng() % 3;
                a = value;
            }
            for (int& f : reps) f = static_cast<int>(rng() % 6) - 1;

            for (int variant = 0; variant < 4; ++variant) {
                if (trial % 2 == 0) {
                    std::fill(p1.begin(), p1.end(), poison);
                    std::fill(p2.begin(), p2.end(), -poison);
                }
                unsigned long expected = 0;
                if (variant >= 2 && strtLen == 0) {
                    expected = 1;
                } else {
                    const int start = variant == 2 ? strtLen : m;
                    for (int j = start; j <= m; ++j) {
                        expected += BruteCount(
                            n, j, allowed, reps, 0, variant != 0
                        );
                    }
                }
                mpz_class result = poison;
                counters[variant](result, p1, p2, n, m, allowed, reps, strtLen);
                const std::string inputs = " trial=" +
                    std::to_string(trial) + " variant=" +
                    std::to_string(variant);
                if (mpz_cmp_ui(result.get_mpz_t(), expected) != 0) {
                    throw std::runtime_error(
                        "Reused buffer count mismatch" + inputs
                    );
                }
                const bool sameStorage =
                    (p1.data() == first && p2.data() == second) ||
                    (p1.data() == second && p2.data() == first);
                if (!sameStorage || p1.size() != size || p2.size() != size) {
                    throw std::runtime_error(
                        "Counting replaced buffer storage" + inputs
                    );
                }
            }
        }

        expect_true(p1.size() == size);
        expect_true(p2.size() == size);
    }

    test_that("ResetP2D clears values without changing row sizes") {
        mpz_class poison = 1;
        poison <<= 128;
        std::vector<std::vector<mpz_class>> p2d = {{poison}, {-poison, poison}};
        ResetP2D(p2d);
        const std::vector<std::vector<mpz_class>> expected = {{0}, {0, 0}};
        // Compare before Catch stringifies operands: streaming mpz_class
        // would introduce a libgmpxx dependency just for diagnostics.
        const bool cleared = p2d == expected;
        expect_true(cleared);
    }
}
