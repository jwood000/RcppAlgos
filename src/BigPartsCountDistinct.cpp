#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountSection.h"
#include <vector>

void CountPartsDistinctLenCap(mpz_class &res, std::vector<mpz_class> &p1,
                              std::vector<mpz_class> &p2, int n, int m,
                              int cap, int strtLen) {

    if (cap > n) cap = n;
    const int limit = (cap * m) - ((m - 1) * m) / 2;

    if (m > n || cap < m) {
        res = 0;
    } else if (m == n) {
        if (n == 1 && cap >= 1) {
            res = 1;
        } else {
            res = 0;
        }
    } else if (m == 1) {
        if (cap >= n) {
            res = 1;
        } else {
            res = 0;
        }
    } else if (limit <= n) {
        if (limit == n) {
            res = 1;
        } else {
            res = 0;
        }
    } else {
        const int width = n + 1;
        const int maxSize = (cap + 1) * width;

        std::fill_n(p1.begin(), maxSize, 0u);

        for (int i = 1; i < width; ++i) {
            for (int j = i; j <= cap; ++j) {
                p1[j * width + i] = 1;
            }
        }

        for (int i = 2; i <= m; ++i) {
            if (i % 2) {
                std::fill_n(p1.begin(), maxSize, 0);

                for (int j = width; j < maxSize; j += width) {
                    for (int k = i, j1 = j - width; k < width; ++k) {
                        p1[j + k] = p2[j1 + k - i] + p1[j1 + k - i];
                    }
                }
            } else {
                std::fill_n(p2.begin(), maxSize, 0);

                for (int j = width; j < maxSize; j += width) {
                    for (int k = i, j1 = j - width; k < width; ++k) {
                        p2[j + k] = p1[j1 + k - i] + p2[j1 + k - i];
                    }
                }
            }
        }

        if (m % 2) {
            res = p1[maxSize - 1];
        } else {
            res = p2[maxSize - 1];
        }
    }
}

void CountPartsDistinctLen(mpz_class &res, std::vector<mpz_class> &p1,
                           std::vector<mpz_class> &p2, int n, int m,
                           int cap, int strtLen) {

    const int max_width = GetMaxWidth(n);

    if (m == 0 && n == 0) {
        res = 1;
    } else if (m == 0) {
        res = 0;
    } else if (m > max_width) {
        res = 0;
    } else if (m < 2) {
        res = 1;
    } else if (m == 2) {
        int n1div2 = (n - 1) / 2;
        res = n1div2;
    } else if (m == 3) {
        mpz_class mpzN(n - 3);
        SumSection(mpzN, res);
    } else {
        const int limit = (m == GetMaxWidth(n + 1)) ? m - 1 : m;

        if (n <= typeSwitchBnd) {
            for (int i = 6; i <= n; ++i) {
                p1[i] = static_cast<double>(SumSection(i - 3));
            }
        } else {
            for (int i = 6; i < typeSwitchBnd; ++i) {
                p1[i] = static_cast<double>(SumSection(i - 3));
            }

            mpz_class tempN;

            for (int i = typeSwitchBnd; i <= n; ++i) {
                tempN = i - 3;
                SumSection(tempN, p1[i]);
            }
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
            res = p2[n - m];
        } else if (m > limit) {
            res = p1[n - m];
        } else if (m % 2) {
            res = p1[n];
        } else {
            res = p2[n];
        }
    }
}

void CountPartsDistinct(mpz_class &res, int n, int m, int cap, int strtLen) {

    std::vector<mpz_class> qq(n + 1);

    qq[0] = 1;
    qq[1] = 1;

    for(int i = 2 ; i <= n; ++i) {
        for (int s = 1, f = 5, r = 2; i >= r; r += f, f += 3, s *= -1) {
            if (s > 0) {
                qq[i] += qq[i - r];

                if (i == r * 2) {
                    --qq[i];
                }
            } else {
                qq[i] -= qq[i - r];

                if (i == r * 2) {
                    ++qq[i];
                }
            }
        }

        for (int s = 1, f = 4, r = 1; i >= r; r += f, f += 3, s *= -1) {
            if (s > 0) {
                qq[i] += qq[i - r];

                if (i == r * 2) {
                    --qq[i];
                }
            } else {
                qq[i] -= qq[i - r];

                if (i == r * 2) {
                    ++qq[i];
                }
            }
        }
    }

    res = qq[n];
}

void CountPartsDistinctMultiZero(mpz_class &res, std::vector<mpz_class> &p1,
                                 std::vector<mpz_class> &p2, int n, int m,
                                 int cap, int strtLen) {

    mpz_class temp;
    res = 0;

    for (int i = strtLen; i <= m; ++i) {
        CountPartsDistinctLen(temp, p1, p2, n, i, cap, strtLen);
        res += temp;
    }
}

void CountPartsDistinctCapMZ(mpz_class &res, std::vector<mpz_class> &p1,
                             std::vector<mpz_class> &p2, int n, int m,
                             int cap, int strtLen) {

    mpz_class temp;
    res = 0;

    for (int i = strtLen; i <= m; ++i) {
        CountPartsDistinctLenCap(temp, p1, p2, n, i, cap, strtLen);
        res += temp;
    }
}
