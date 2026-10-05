#include "Combinations/BigComboCount.h"
#include <algorithm>
#include <gmpxx.h>
#include <vector>

// See notes in PartitionsCountMultiset.cpp for general invariants.
//
// The algorithms below differ slightly from their `double` counterparts because
// of optimizations specific to `mpz_class`. The most noticeable difference is
// how the working vectors are reset between iterations. In the `double` case,
// creating a fresh vector each iteration is relatively cheap. In the GMP
// partition algorithm, we instead make judicious use of sliding windows and
// buffer swapping. The GMP composition algorithm cannot use the same sliding
// window recurrence, so it more closely resembles the `double` implementation,
// aside from the use of persistent buffers and swapping.
//
// Flattened (m + 1) x (n + 1):
// dp[j][s] stored at dp[j * (n + 1) + s]
static inline std::size_t idx_b(std::size_t j, std::size_t s, std::size_t n1) {
    return j * n1 + s;
}

void CountPartsMultisetWorker(
    mpz_class &result, int n, int m,
    const std::vector<int>& allowed,
    const std::vector<int>& Reps,
    int strtLen = 0
) {

    const int k = static_cast<int>(allowed.size());
    const std::size_t n1 = static_cast<std::size_t>(n) + 1;
    const std::size_t size = (static_cast<std::size_t>(m) + 1) * n1;

    std::vector<mpz_class> curr(size, 0);
    std::vector<mpz_class> next(size, 0);

    curr[idx_b(0, 0, n1)] = 1;
    int maxParts = 0;

    // allowed is sorted and positionally aligned with Reps.
    // Ignore non-positive entries; once a > n, no later value can contribute.
    for (int i = 0; i < k; ++i) {
        const int a = allowed[i];
        if (a < 1) continue;
        if (a > n) break;

        const int f = std::min({Reps[i], m, n / a});
        if (f < 1) continue;

        const int window = f + 1;

        if (f >= m - maxParts) {
            maxParts = m;
        } else {
            maxParts += f;
        }

        for (int j = 0; j <= maxParts; ++j) {
            const int s_max = std::min(n, j * a);
            const std::size_t row = j * n1;

            for (int s = j; s <= s_max; ++s) {
                const std::size_t pos = row + s;

                // t = 0
                next[pos] = curr[pos];

                // Add the sliding window from the previous j.
                if (j > 0 && s >= a) {
                    next[pos] += next[idx_b(j - 1, s - a, n1)];
                }

                // Remove the term that has fallen outside
                // the allowed multiplicity window.
                if (j >= window && s >= window * a) {
                    next[pos] -= curr[idx_b(j - window, s - window * a, n1)];
                }
            }
        }

        curr.swap(next);
    }

    if (strtLen) {
        result = 0;

        for (int j = strtLen; j <= m; ++j) {
            result += curr[idx_b(j, n, n1)];
        }
    } else {
        result = curr[idx_b(m, n, n1)];
    }
}

void CountCompsMultisetWorker(
    mpz_class &result, int n, int m,
    const std::vector<int>& allowed,
    const std::vector<int>& Reps,
    int strtLen = 0
) {

    const int k = static_cast<int>(allowed.size());
    const std::size_t n1 = static_cast<std::size_t>(n) + 1;
    const std::size_t size = (static_cast<std::size_t>(m) + 1) * n1;

    std::vector<mpz_class> curr(size, 0);
    std::vector<mpz_class> next(size, 0);

    curr[idx_b(0, 0, n1)] = 1;
    int maxParts = 0;
    mpz_class mult = 1;

    // allowed is sorted and positionally aligned with Reps.
    // Ignore non-positive entries; once a > n, no later value can contribute.
    for (int i = 0; i < k; ++i) {
        const int a = allowed[i];
        if (a < 1) continue;
        if (a > n) break;

        const int f = std::min({Reps[i], m, n / a});
        if (f < 1) continue;

        const int oldMaxParts = maxParts;

        if (f >= m - maxParts) {
            maxParts = m;
        } else {
            maxParts += f;
        }

        const int t_max = std::min(f, m);
        std::fill(next.begin(), next.end(), 0);

        for (int t = 0; t <= t_max; ++t) {
            const int shift = t * a;
            if (shift > n) break;

            const int j_max = std::min(oldMaxParts, m - t);

            // dp[j + t][s + shift] += dp_old[j][s]
            for (int j = 0; j <= j_max; ++j) {
                const std::size_t row_src = idx_b(j, 0, n1);
                const std::size_t row_dst = idx_b(j + t, 0, n1);
                const int s_max = std::min(n - shift, j * (a - 1));
                nChooseKGmp(mult, j + t, t);

                for (int s = j; s <= s_max; ++s) {
                    // Avoid mpz_class temporaries in this hot loop.
                    if (cmp(curr[row_src + s], 0) > 0) {
                        mpz_addmul(
                            next[row_dst + s + shift].get_mpz_t(),
                            curr[row_src + s].get_mpz_t(),
                            mult.get_mpz_t()
                        );
                    }
                }
            }
        }

        curr.swap(next);
    }

    if (strtLen) {
        result = 0;

        for (int j = strtLen; j <= m; ++j) {
            result += curr[idx_b(j, n, n1)];
        }
    } else {
        result = curr[idx_b(m, n, n1)];
    }
}

void CountPartsMultiset(
    mpz_class &res, int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
) {
    return CountPartsMultisetWorker(res, n, m, allowed, Reps);
}

void CountCompsMultiset(
    mpz_class &res, int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
) {
    return CountCompsMultisetWorker(res, n, m, allowed, Reps);
}

void CountCompsMultisetZero(
    mpz_class &res, int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        res = 1;
    } else {
        CountCompsMultisetWorker(res, n, m, allowed, Reps, strtLen);
    }
}

void CountCompsMultisetWeak(
    mpz_class &res, int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        res = 1;
    } else {
        return CountCompsMultisetWorker(res, n, m, allowed, Reps);
    }
}
