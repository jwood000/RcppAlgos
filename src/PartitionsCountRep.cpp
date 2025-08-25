#include "Partitions/PartitionsCountSection.h"
#include "Combinations/ComboCount.h"
#include <algorithm>  // std::fill; std::min
#include <cmath>      // std::floor

double CountPartsRepLenRstrctd(
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    std::vector<std::vector<double>> p(m + 1, std::vector<double>(n + 1, 0));
    p[0][0] = 1; // one way to make 0 with 0 parts

    for (int num : allowed) {
        for (int j = 1; j <= m; ++j) {
            for (int s = num; s <= n; ++s) {
                p[j][s] += p[j - 1][s - num];
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
