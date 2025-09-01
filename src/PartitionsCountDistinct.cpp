#include "Partitions/PartitionsCountSection.h"
#include "Permutations/PermuteCount.h"
#include <algorithm>  // std::count_if
#include <numeric>    // std::iota

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
