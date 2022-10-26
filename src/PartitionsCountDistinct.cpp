#include "Partitions/PartitionsCountSection.h"
#include "Permutations/PermuteCount.h"
#include <algorithm>  // std::count_if
#include <numeric>
#include <cmath>

double CountPartsDistinctLenCap(int n, int m, int cap, int strtLen) {

    if (cap > n) cap = n;
    if (m > n || cap < m) return 0;

    if (m == n) {
        if (n == 1 && cap >= 1) {
            return 1;
        } else {
            return 0;
        }
    }

    if (m == 1) {
        if (cap >= n) {
            return 1;
        } else {
            return 0;
        }
    }

    // Ensure max is large enough given the width
    //
    // Below, we have an expression that represents the
    // absolute maximum value we could obtain with a
    // given max and width:
    //
    // max + (max - 1) + (max - 2) + ... + (max - (m - 1))
    //
    // (max * m) - (0 + 1 + 2 + ... + (m - 1))
    //
    // (max * m) - ((m - 1) * m) / 2

    CheckMultIsInt(cap, m);
    CheckMultIsInt(m - 1, m);
    const int limit = (cap * m) - ((m - 1) * m) / 2;

    if (limit <= n) {
        if (limit == n) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    const int width = n + 1;
    CheckMultIsInt(cap + 1, width);
    const int maxSize = (cap + 1) * width;

    std::vector<double> p1(maxSize);
    std::vector<double> p2(maxSize);

    for (int i = 1; i < width; ++i) {
        for (int j = i; j <= cap; ++j) {
            p1[j * width + i] = 1;
        }
    }

    for (int i = 2; i <= m; ++i) {
        if (i % 2) {
            std::fill(p1.begin(), p1.end(), 0);

            for (int j = width; j < maxSize; j += width) {
                for (int k = i, j1 = j - width; k < width; ++k) {
                    p1[j + k] = p2[j1 + k - i] + p1[j1 + k - i];
                }
            }
        } else {
            std::fill(p2.begin(), p2.end(), 0);

            for (int j = width; j < maxSize; j += width) {
                for (int k = i, j1 = j - width; k < width; ++k) {
                    p2[j + k] = p1[j1 + k - i] + p2[j1 + k - i];
                }
            }
        }
    }

    return (m % 2) ? p1.back() : p2.back();
}

double CountPartsDistinctLen(int n, int m, int cap, int strtLen) {

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
    } else {
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
}

// Credit to Robin K. S. Hankin, author of the excellent partitions package.
// From the partitions.c, here are Hankin's comments for c_numbdiffparts:
//      "the recursion on p826 of Abramowitz and Stegun"
double CountPartsDistinct(int n, int m, int cap, int strtLen) {

    std::vector<double> qq(n + 1);
    qq[0] = 1;
    qq[1] = 1;

    for(int i = 2 ; i <= n; ++i) {
        for (int s = 1, f = 5, r = 2; i >= r; r += f, f += 3, s *= -1) {
            qq[i] += s * qq[i - r];
            if(i == r * 2) {qq[i] -= s;}
        }

        for (int s = 1, f = 4, r = 1; i >= r; r += f, f += 3, s *= -1) {
            qq[i] += s * qq[i - r];
            if(i == r * 2) {qq[i] -= s;}
        }
    }

    return qq.back();
}

double CountPartsDistinctMultiZero(int n, int m, int cap, int strtLen) {

    double count = 0;

    for (int i = strtLen; i <= m; ++i) {
        count += CountPartsDistinctLen(n, i, cap, strtLen);
    }

    return count;
}

double CountPartsDistinctCapMZ(int n, int m, int cap, int strtLen) {

    double count = 0;

    for (int i = strtLen; i <= m; ++i) {
        count += CountPartsDistinctLenCap(n, i, cap, strtLen);
    }

    return count;
}

double CountPartsPermDistinct(const std::vector<int> &z,
                              int tar, int width, bool includeZero) {

    double res = 0;

    if (includeZero) {
        const int strtLen = std::count_if(z.cbegin(), z.cend(),
                                          [](int i){return i > 0;});

        if (strtLen == 0) {
            // This means that z contains only zeros
            res = 1;
        } else {
            std::vector<int> count(width);
            std::iota(count.begin(), count.begin() + strtLen, 1);

            for (int i = strtLen; i <= width; ++i) {
                count[i - 1] = i;
                res += (CountPartsDistinctLen(tar, i, tar, tar) *
                        NumPermsWithRep(count));
            }
        }
    } else {
        res = CountPartsDistinctLen(tar, width, tar, tar) *
              NumPermsNoRep(width, width);
    }

    return res;
}

double CountPartsPermDistinctCap(const std::vector<int> &z, int cap,
                                 int tar, int width, bool includeZero) {

    double res = 0;

    if (includeZero) {
        const int strtLen = std::count_if(z.cbegin(), z.cend(),
                                           [](int i){return i > 0;});
        if (strtLen == 0) {
            // This means that z contains only zeros
            res = 1;
        } else {
            std::vector<int> count(width);
            std::iota(count.begin(), count.begin() + strtLen, 1);

            for (int i = strtLen; i <= width; ++i) {
                count[i - 1] = i;
                res += (CountPartsDistinctLenCap(tar, i, cap, tar) *
                        NumPermsWithRep(count));
            }
        }
    } else {
        res = CountPartsDistinctLenCap(tar, width, cap, tar) *
              NumPermsNoRep(width, width);
    }

    return res;
}
