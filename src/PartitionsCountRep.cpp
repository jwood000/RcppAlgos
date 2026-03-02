#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountRep.h"
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

// Returns the count as a double. Internally, we compute the result using the
// GMP (mpz_class) routine and convert at the end.
//
// Rationale:
// This count is based on an inclusion–exclusion formula. Individual terms
// (e.g. choose(n - 1, m - 1)) can exceed 2^53 - 1 even when the final result
// is smaller due to cancellation. If we performed the inclusion–exclusion
// arithmetic directly in double, those large intermediate values would lose
// integer precision and the final result could be off by one or more.
//
// By computing the sum exactly with mpz_class first and only converting the
// final value to double, we ensure correctness whenever the true result is
// representable exactly in double (<= 2^53 - 1). This preserves the double
// interface while avoiding floating-point drift in razor-edge cases.
double CountCompsRepLenCap(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    mpz_class res;
    CountCompsRepLenCap(res, n, m, allowed);
    return res.get_d();
}

double CountCompsRepCapZNotWk(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    mpz_class res;
    CountCompsRepCapZNotWk(res, n, m, allowed);
    return res.get_d();
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
