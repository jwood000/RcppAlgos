#include "Partitions/MultisetCountEngine.h"
#include "Combinations/ComboCount.h"
#include <algorithm>

// Multiset counting invariants:
//
// * `allowed` contains the mapped positive part values used by the DP.
//   It is sorted in strictly increasing order and is positionally aligned
//   with `Reps`, so `Reps[i]` is the multiplicity of `allowed[i]`.
//
// * `allowed.size()` will always be less than or equal to `Reps.size()`. When
//   the partition case is mapped, the maximum value, held in `part.cap`, may
//   be less than `max(v)`. For this reason, `k` is determined by `allowed`.
//
// * Zero-valued inputs are normalized before this stage. `allowed` refers to
//   mapped DP values, not the user's original values.
//
// * Cases in which all multiplicities are one are routed to the distinct
//   counting path before reaching these workers.
//
// * `strtLen == 0` is permitted because ranking may request a count from a
//   state with no positive parts remaining. The caller is responsible for
//   establishing that the state is valid; this routine does not re-validate it.

// dp[j][s] = # ways to write s using exactly j parts,
//            where part size allowed[i] can be used at most Reps[i] times.
// For compositions, each transition is weighted by the number of
// distinct interleavings of the newly added equal parts.

template <bool IsComp>
double MultisetCountEngine<IsComp>::Count(
    int n, int m,
    const std::vector<int>& allowed,
    const std::vector<int>& Reps,
    int strtLen
) {

    const int k = static_cast<int>(allowed.size());
    const std::size_t n1   = static_cast<std::size_t>(n) + 1;
    const std::size_t size = (static_cast<std::size_t>(m) + 1) * n1;

    dp.assign(size, 0.0);
    dp_old.resize(size);

    dp[idx(0, 0, n1)] = 1.0; // 0 parts sum to 0 in one way
    int maxParts = 0;

    // allowed is sorted and positionally aligned with Reps.
    // Ignore non-positive entries; once a > n, no later value can contribute.
    for (int i = 0; i < k; ++i) {
        const int a = allowed[i];
        if (a < 1) continue;
        if (a > n) break;

        const int f = std::min({Reps[i], m, n / a});
        if (f < 1) continue;

        // Snapshot before using part size a
        std::copy(dp.begin(), dp.begin() + size, dp_old.begin());
        const int t_max = std::min(f, m);
        const int oldMaxParts = maxParts;

        if (f >= m - maxParts) {
            maxParts = m;
        } else {
            maxParts += f;
        }

        for (int t = 1; t <= t_max; ++t) {
            const int shift = t * a;
            if (shift > n) break;

            const int j_max = std::min(oldMaxParts, m - t);

            // dp[j + t][s + shift] += dp_old[j][s]
            for (int j = 0; j <= j_max; ++j) {
                const std::size_t row_src = idx(j, 0, n1);
                const std::size_t row_dst = idx(j + t, 0, n1);
                const double mult = IsComp ? nChooseK(j + t, t) : 1.0;

                for (int s = 0; s <= n - shift; ++s) {
                    const double src = dp_old[row_src + s];

                    if (src != 0.0) {
                        dp[row_dst + s + shift] += src * mult;
                    }
                }
            }
        }
    }

    if (strtLen) {
        double result = 0.0;

        for (int j = strtLen; j <= m; ++j) {
            result += dp[idx(j, n, n1)];
        }

        return result;
    } else {
        return dp[idx(m, n, n1)];
    }
}

template double MultisetCountEngine<true>::Count(
    int, int, const std::vector<int>&, const std::vector<int>&, int
);

template double MultisetCountEngine<false>::Count(
    int, int, const std::vector<int>&, const std::vector<int>&, int
);
