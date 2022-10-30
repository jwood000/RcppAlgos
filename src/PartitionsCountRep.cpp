#include "Partitions/PartitionsCountSection.h"
#include "Combinations/ComboCount.h"
#include <algorithm>  // std::fill; std::min
#include <cmath>      // std::floor

double CountPartsRepLenCap(int n, int m, int cap, int strtLen) {

    if (cap > n) cap = n;
    CheckMultIsInt(cap, m);
    if (cap * m < n || n < m) return 0;
    if (cap * m == n || n <= m + 1) return 1;
    if (m < 2) return m;

    if (m == 2) {
        CheckMultIsInt(2, cap);
        if (cap * 2 >= n) {
            cap = std::min(cap, n - 1);
            return n / m - (n - 1 - cap);
        } else {
            return 0;
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
                    p1[j + k] = p2[j + k - 1] + p1[j1 + k - i];
                }
            }
        } else {
            std::fill(p2.begin(), p2.end(), 0);

            for (int j = width; j < maxSize; j += width) {
                for (int k = i, j1 = j - width; k < width; ++k) {
                    p2[j + k] = p1[j + k - 1] + p2[j1 + k - i];
                }
            }
        }
    }

    return (m % 2) ? p1.back() : p2.back();
}

double CountPartsRepLen(int n, int m, int cap, int strtLen) {

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
    } else {
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
}

// Similar to CountPartsDistinct
double CountPartsRep(int n, int m, int cap, int strtLen) {

    if (n < 2) {
        return 1.0;
    }

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

double CountCompsRepLen(int n, int m, int cap, int strtLen) {
    return nChooseK(n - 1, m - 1);
}

double CountCompsRepZero(int n, int m, int cap, int strtLen) {

    if (n == m) return std::pow(2.0, static_cast<double>(n - 1));

    double res = 0;

    for (int i = 1; i <= m; ++i) {
        res += nChooseK(n - 1, i - 1);
    }

    return res;
}
