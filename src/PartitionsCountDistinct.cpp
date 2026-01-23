#include "Partitions/PartitionsCountSection.h"
#include "Permutations/PermuteCount.h"
#include <algorithm>  // std::count_if
#include <numeric>    // std::iota

// UpdateAllowed
// -------------
// Maintains two coupled data structures used by nthCompsDistinct:
//
//   1) mask[v] = 1 if value v is already used in the partial solution,
//                 0 if it is still available.
//   2) allowed  = compact list of still-available candidate values for the
//                  remaining positions, after accounting for a feasibility
//                  bound.
//
// Invariants assumed by nthCompsDistinct:
//   - Parts are distinct, positive integers.
//   - We are unranking in lexicographic order.
//   - At position i we just fixed a value new_val (and are replacing the
//     previous "current" value cur_val in the bookkeeping).
//   - The remaining suffix has (width - i - 1) positions.
//   - We only keep values in allowed that are feasible given the remaining sum.
//
// What this function computes:
//   - It updates the mask to mark new_val as used and cur_val as unused.
//   - It computes a conservative upper bound last_val for the next candidates.
//     Any candidate > last_val cannot work because even choosing the smallest
//     unused values for the remaining slots would overshoot the sum constraint.
//
// The partial_sum parameter is the running sum of the already-fixed prefix
// including the current position after it is updated by the caller.
//
void UpdateAllowed(
    std::vector<char> &mask, std::vector<int> &allowed, int i,
    int new_val, int width, int n, int cur_val, int partial_sum
) {

    // "Move" the used marker from cur_val -> new_val.
    // This allows the caller to sweep through potential candidates and
    // cheaply update feasibility.
    mask[cur_val] = 0;
    mask[new_val] = 1;

    // Add the minimum possible contribution of the remaining suffix
    // (excluding the final slot, which will be forced at the very end).
    //
    // For each remaining internal position, we greedily pick the smallest
    // unused positive integer. If even this minimal tail makes the total
    // exceed n, then no completion is possible.
    //
    // Loop structure:
    //   - j indexes future positions
    //   - k is the current candidate value, advanced until unused.
    for (int j = i + 1, k = 1; j < (width - 1); ++j, ++k) {
        while (mask[k]) {
            ++k;
        }

        partial_sum += k;
    }

    // Given the minimal tail contribution computed above, the next selected
    // value cannot exceed (n - partial_sum). Also, it cannot exceed cap
    // (mask.size() - 1).
    int j = 0;
    const int last_val = std::min(
        n - partial_sum,
        static_cast<int>(mask.size()) - 1
    );

    // Build the compact allowed list = all unused values in [1..last_val].
    // (Anything larger is unfeasible by construction.)
    for (int v = 1; v <= last_val; ++v) {
        if (!mask[v]) {
            allowed[j] = v;
            ++j;
        }
    }

    // Zero out the remainder as padding (so loops can skip if (num)).
    std::fill(allowed.begin() + j, allowed.end(), 0);
}

// CountPartsDistLenRstrctd
// ------------------------
// Counts the number of ways to choose exactly m distinct values from allowed
// such that they sum to n.
//
// Interpretation:
//   - This is a subset-counting DP (order does NOT matter).
//   - Each element of allowed may be used at most once.
//   - allowed may contain zeros as padding; zeros are ignored (if (num)).
//
// DP definition:
//   dp[j][s] = number of ways to pick j distinct values that sum to s
//
// Transition (classic 0/1 knapsack style, iterating backwards to avoid reuse):
//   dp[j][s] += dp[j - 1][s - num]
//
// Complexity:
//   O(|allowed|  m  n) time, O(m * n) memory.
//
// Through extensive testing, this 2D approach is much more efficient.
// This is probably due to allocating a huge, mostly empty and unused,
// vector for the 3D case.
double CountPartsDistLenRstrctd(
    int n, int m, const std::vector<int> &allowed, int strtLen = 0
) {

    std::vector<std::vector<double>> dp(m + 1, std::vector<double>(n + 1, 0));
    dp[0][0] = 1; // one way to make 0 with 0 parts

    for (int num: allowed) {
        if (num) {
            for (int j = m; j >= 1; --j) {
                for (int s = n; s >= num; --s) {
                    dp[j][s] += dp[j - 1][s - num];
                }
            }
        }
    }

    return dp[m][n];
}

double CountPartsDistinctLen(
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
) {

    const int max_width = GetMaxWidth(n);

    if (m == 0) {
        return (n == 0) ? 1.0 : 0.0;
    } else if (m > max_width) {
        return 0.0;
    } else if (m < 2) {
        return 1.0;
    } else if (m == 2) {
        return (n - 1) / 2;
    } else if (m == 3) {
        const double res = SumSection(n - 3);
        return(res);
    }

    const int limit = (m == GetMaxWidth(n + 1)) ? m - 1 : m;
    std::vector<double> p1(n + 1);
    std::vector<double> p2(n + 1);

    for (int i = 6; i <= n; ++i) {
        p1[i] = SumSection(i - 3);
    }

    for (int i = 4; i <= limit; ++i) {
        const int m1 = ((i + 1) * i) / 2;
        const int m2 = m1 + i;

        if (i % 2) {
            for (int j = m1; j < m2; ++j) {
                p1[j] = p2[j - i];
            }

            for (int j = m2; j <= n; ++j) {
                p1[j] = p2[j - i] + p1[j - i];
            }
        } else {
            for (int j = m1; j < m2; ++j) {
                p2[j] = p1[j - i];
            }

            for (int j = m2; j <= n; ++j) {
                p2[j] = p2[j - i] + p1[j - i];
            }
        }
    }

    if (m > limit && m % 2) {
        return p2[n - m];
    } else if (m > limit) {
        return p1[n - m];
    } else if (m % 2) {
        return p1.back();
    } else {
        return p2.back();
    }
}

// Credit to Robin K. S. Hankin, author of the excellent partitions package.
// From the partitions.c, here are Hankin's comments for c_numbdiffparts:
//      "the recursion on p826 of Abramowitz and Stegun"
double CountPartsDistinct(
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
) {

    std::vector<double> qq(n + 1);
    qq[0] = 1;
    qq[1] = 1;

    for(int i = 2 ; i <= n; ++i) {
        for (int s = 1, f = 5, r = 2; i >= r; r += f, f += 3, s *= -1) {
            qq[i] += s * qq[i - r];

            if (i == r * 2) {
                qq[i] -= s;
            }
        }

        for (int s = 1, f = 4, r = 1; i >= r; r += f, f += 3, s *= -1) {
            qq[i] += s * qq[i - r];

            if (i == r * 2) {
                qq[i] -= s;
            }
        }
    }

    return qq.back();
}

double CountPartsDistinctMultiZero(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    double count = 0;

    for (int i = strtLen; i <= m; ++i) {
        count += CountPartsDistinctLen(n, i);
    }

    return count;
}

double CountPartsDistinctRstrctdMZ(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    double count = 0;

    for (int i = strtLen; i <= m; ++i) {
        count += CountPartsDistLenRstrctd(n, i, allowed);
    }

    return count;
}

// CountCompDistLenRstrctd
// -----------------------
// Counts the number of compositions (order matters) of length m whose parts
// are distinct values chosen from allowed and sum to n.
//
// Key idea:
//   - First count distinct subsets of size m summing to n (orderless) via
//     CountPartsDistLenRstrctd.
//   - Then convert each subset to all possible orderings.
//   - Because all chosen values are distinct, the number of orderings is m!.
//     (NumPermsNoRep(m, m) is assumed to be m!.)
//
double CountCompDistLenRstrctd(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsDistLenRstrctd(n, m, allowed) * NumPermsNoRep(m, m);
}

double CountPartsPermDistinctRstrctdMZ(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        return 1;
    } else {
        double res = 0;
        std::vector<int> v(m);
        std::iota(v.begin(), v.begin() + strtLen, 1);

        for (int i = strtLen; i <= m; ++i) {
            v[i - 1] = i;
            res += (CountPartsDistLenRstrctd(n, i, allowed) *
                NumPermsWithRep(v));
        }

        return res;
    }
}

double CountCompsDistinctLen(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return CountPartsDistinctLen(n, m) * NumPermsNoRep(m, m);
}

double CountCompsDistinctMultiZero(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        return 1;
    } else {
        double res = 0;
        double nPerm = NumPermsNoRep(strtLen, strtLen);

        for (int i = strtLen; i <= m; ++i) {
            res += CountPartsDistinctLen(n, i) * nPerm;
            nPerm *= (i + 1);
        }

        return res;
    }
}

double CountCompsDistinctMZWeak(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        return 1;
    } else {
        double res = 0;
        double nPerm = NumPermsNoRep(m, m) /
            NumPermsNoRep(m - strtLen, m - strtLen);

        for (int i = strtLen; i <= m; ++i) {
            res += CountPartsDistinctLen(n, i) * nPerm;
            nPerm *= (m - i);
        }

        return res;
    }
}

double CountCompsDistinctRstrctdMZ(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        return 1;
    } else {
        double res = 0;
        double nPerm = NumPermsNoRep(strtLen, strtLen);

        for (int i = strtLen; i <= m; ++i) {
            res += CountPartsDistLenRstrctd(n, i, allowed) * nPerm;
            nPerm *= (i + 1);
        }

        return res;
    }
}

double CountCompsDistinctRstrctdMZWeak(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        return 1;
    } else {
        double res = 0;
        double nPerm = NumPermsNoRep(m, m) /
            NumPermsNoRep(m - strtLen, m - strtLen);

        for (int i = strtLen; i <= m; ++i) {
            res += CountPartsDistLenRstrctd(n, i, allowed) * nPerm;
            nPerm *= (m - i);
        }

        return res;
    }
}
