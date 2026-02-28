#include "Partitions/PartitionsCountSection.h"
#include "Combinations/ComboCount.h"
#include <algorithm>  // std::fill; std::min
#include <cmath>      // std::floor

// CountPartsRepLenRstrctd
//
// PRECONDITION: `allowed` encodes coin types as strictly increasing positive
// values, optionally followed by trailing zeros as padding. Zeros are ignored.
//
// Counts the number of partitions of length `m` that sum to `n`
// where:
//
//   • Parts are drawn from the set `allowed`.
//   • Repetition of values IS allowed.
//   • Order does NOT matter (i.e. partitions, not compositions).
//
// IMPORTANT: Semantics of `allowed`
// ----------------------------------
// `allowed` represents the set of admissible *values*, not a multiset
// of consumable elements. Each value in `allowed` may be used multiple
// times in the partition (subject only to the length constraint `m`).
//
// This is therefore the classic fixed-length coin-change / combinations-
// with-repetition problem.
//
// DP Recurrence
// -------------
// Let dp[j][s] denote the number of ways to form sum `s` using exactly
// `j` parts from `allowed`, with repetition allowed and order ignored.
//
// The recurrence:
//
//   dp[j][s] += dp[j - 1][s - num]
//
// for each `num ∈ allowed`
//
// is correct because:
//
//   • We build partitions by appending `num` as the final part.
//   • dp[j - 1][s - num] already counts all partitions of length j-1
//     summing to s-num.
//   • Since order does not matter and we iterate values consistently,
//     each unordered partition is counted exactly once.
//   • Repetition is allowed because the same `num` may be reused in
//     subsequent transitions; values are not “consumed”.
//
// This is NOT a 0/1 subset-sum recurrence. Values are reusable.
//
// Validation
// ----------
// Verified against brute-force enumeration:
//
//   comboGeneral(cap, m, repetition = TRUE, constraintFun = "sum")
//
// bucketed by sum, and compared to:
//
//   partitionsCount(cap, m, repetition = TRUE, target = s)
//
// across all feasible target sums. Results match exactly.
//
double CountPartsRepLenRstrctd(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    CheckAllowedInvariant(allowed);
    std::vector<std::vector<double>> p(m + 1, std::vector<double>(n + 1, 0));
    p[0][0] = 1; // one way to make 0 with 0 parts

    for (int num : allowed) {
        if (num) {
            for (int j = 1; j <= m; ++j) {
                for (int s = num; s <= n; ++s) {
                    p[j][s] += p[j - 1][s - num];
                }
            }
        }
    }

    return p[m][n];
}

double CountPartsRepLen(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (m == 0) {
        return (n == 0) ? 1.0 : 0.0;
    } else if (n < m) {
        return 0.0;
    } else if (n == m) {
        return 1.0;
    } else if (m < 2) {
        return 1.0;
    } else if (n - m == 1) {
        return 1.0;
    } else if (m == 2) {
        return std::floor(n / 2);
    }

    // If n > 3, we have the following:
    // 1st part: n - 3 1's followed by a 3
    // 2nd part: n - 4 1's followed by two 2's
    //
    // N.B. We have taken care of every case where n <= 2 above
    // i.e. n = 2, m = 2; n = 2, m = 1; n = 1, m = 1; n = 0, m = 0
    if (n - m == 2) {
        return 2.0;
    }

    if (m == 3) {
        const double res = SumSection(n);
        return(res);
    }

    const int limit = std::min(n - m, m);
    CheckMultIsInt(2, m);
    CheckMultIsInt(2, limit);
    n = (n < 2 * m) ? 2 * limit : n;

    std::vector<double> p1(n + 1);
    std::vector<double> p2(n + 1);

    for (int i = 3; i <= n; ++i) {
        p1[i] = SumSection(i);
    }

    for (int i = 4; i <= limit; ++i) {
        const int m2 = i * 2;

        if (i % 2) {
            p1[i] = 1;

            for (int j = i + 1; j < m2; ++j) {
                p1[j] = p2[j - 1];
            }

            for (int j = m2; j <= n; ++j) {
                p1[j] = p2[j - 1] + p1[j - i];
            }
        } else {
            p2[i] = 1;

            for (int j = i + 1; j < m2; ++j) {
                p2[j] = p1[j - 1];
            }

            for (int j = m2; j <= n; ++j) {
                p2[j] = p1[j - 1] + p2[j - i];
            }
        }
    }

    return (limit % 2) ? p1.back() : p2.back();
}

// Similar to CountPartsDistinct
double CountPartsRep(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (n < 2) return 1.0;
    std::vector<double> qq(n + 1);

    qq[0] = 1;
    qq[1] = 1;

    for(int i = 2; i <= n; ++i) {
        for (int s = 1, f = 1, r = 1; i >= r; f += 3, r += f, s *= -1) {
            qq[i] += s * qq[i - r];
        }

        for (int s = 1, f = 2, r = 2; i >= r; f += 3, r += f, s *= -1) {
            qq[i] += s * qq[i - r];
        }
    }

    return qq.back();
}

double CountCompsRepLen(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {
    return nChooseK(n - 1, m - 1);
}

// ******* CountCompsRepLenCap Details ******
//
// This uses the inclusion–exclusion principle to count the number of positive
// compositions of length m that sum to n, with each part bounded above by cap.
//
// First, count all positive compositions without the cap:
//
//     x_1 + x_2 + ... + x_m = n,    x_i >= 1
//
// This is:
//
//     choose(n - 1, m - 1)
//
// Next, subtract compositions where at least one part exceeds the cap.
// Suppose exactly one part violates the cap. Without loss of generality:
//
//     x_1 >= cap + 1
//
// Write:
//
//     x_1 = cap + z_1
//     x_i = z_i   for i >= 2
//
// Then:
//
//     z_1 + z_2 + ... + z_m = n - cap
//
// The number of such compositions is:
//
//     choose(n - cap - 1, m - 1)
//
// and there are choose(m, 1) choices for which part violates the cap.
//
// However, this subtracts too much, because compositions where two parts
// exceed the cap were subtracted twice. To correct this, add back the number
// of compositions where two parts exceed the cap:
//
//     z_1 + ... + z_m = n - 2*cap
//
// giving:
//
//     choose(n - 2*cap - 1, m - 1)
//
// and choose(m, 2) ways to select the violating parts.
//
// Continuing in this way, alternating subtraction and addition, yields:
//
//     sum_{i=0}^{floor((n - m)/cap)}
//         (-1)^i * choose(m, i) * choose(n - i*cap - 1, m - 1)
//
// The upper limit floor((n - m)/cap) arises because each violating part must
// contribute at least cap extra beyond the minimum value of 1.
//
// Count capped positive compositions with repetition.
// PRECONDITIONS (enforced by caller / dispatch):
//   - n >= m >= 1
//   - cap >= 2
//   - n <= cap * m
//   - allowed is a single-element vector containing cap (allowed.size() == 1),
//     or the set [1..cap] defined in PartitionsCount
//
// NOTE: This function does NOT support arbitrary allowed-sets. It counts
//       parts in the standard range [1..cap].
double CountCompsRepLenCap(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    static_cast<void>(strtLen); // intentionally unused

    const int cap = *std::max_element(allowed.cbegin(), allowed.cend());
    const int maxViolations = (n - m) >= 0 ? (n - m) / cap : -1;
    double res = nChooseK(n - 1, m - 1);

    for (int i = 1, sign = -1; i <= maxViolations; ++i, sign *= -1) {
        res += (sign * nChooseK(m, i) * nChooseK(n - i * cap - 1, m - 1));
    }

    return res;
}

// The "Z" means that zero is in the output but not considered.
//
// For weak compositions, in SetPartitionDesign in PartitionsUtils.cpp, if we
// are dealing with weak compositions, we simply map 0 to 1. Of course we must
// adjust the target to reflect this. Counting now is trivial as we don't have
// to treat zero special... we simply take advantage of the mapping:
//
//                                      j zeros
// Original first partition: part = c(0, 0, ..., 0, x_(j + 1), ..., x_n)
//
// sum(part) = x_(j + 1) + ... + x_n = P
//
// New mapping: part' = c(1, 1, ..., 1, x_(j + 1) + 1, ..., x_n + 1)
//
//                j 1's          (n - j) terms incremented by 1
// sum(part') = 1 + 1 + ... + 1 + x_(j + 1) + 1 + ... + x_n + 1 =
//
// Since we have incremented each term we have:
//
// sum(part') = n + P

double CountCompsRepZNotWk(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (n == m) return std::pow(2.0, static_cast<double>(n - 1));

    double res = 0;

    for (int i = 1; i <= m; ++i) {
        res += nChooseK(n - 1, i - 1);
    }

    return res;
}
