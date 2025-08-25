#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountSection.h"
#include "Permutations/BigPermuteCount.h"
#include <algorithm>  // std::count_if
#include <numeric>    // std::iota

void CountPartsDistLenRstrctd(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen = 0
) {

    p2d[0][0] = 1; // one way to make 0 with 0 parts

    for (int num: allowed) {
        for (int j = m; j >= 1; --j) {
            for (int s = n; s >= num; --s) {
                p2d[j][s] += p2d[j - 1][s - num];
            }
        }
    }

    res = p2d[m][n];
}

void CountPartsDistinctLen(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
) {

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

void CountPartsDistinct(
    mpz_class &res, int n, int m, const std::vector<int> &allowed,
    int strtLen
) {

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

void CountPartsDistinctMultiZero(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    mpz_class temp;
    res = 0;

    for (int i = strtLen; i <= m; ++i) {
        CountPartsDistinctLen(temp, p1, p2, n, i);
        res += temp;
    }
}

void CountPartsDistinctRstrctdMZ(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    mpz_class temp;
    res = 0;

    for (int i = strtLen; i <= m; ++i) {
        CountPartsDistLenRstrctd(temp, p2d, n, i, allowed);
        res += temp;
        ResetP2D(p2d);
    }
}

void CountPartsPermDistinctRstrctd(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    mpz_class partsCnt = 0;
    mpz_class permsCnt = 0;

    NumPermsNoRepGmp(permsCnt, m, m);
    CountPartsDistLenRstrctd(partsCnt, p2d, n, m, allowed);

    res = (partsCnt * permsCnt);
}

void CountPartsPermDistinctRstrctdMZ(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        res = 1;
    } else {
        res = 0;
        mpz_class partsCnt = 0;
        mpz_class permsCnt = 0;
        NumPermsNoRepGmp(permsCnt, strtLen, strtLen);

        for (int i = strtLen; i <= m; ++i) {
            CountPartsDistLenRstrctd(partsCnt, p2d, n, i, allowed);
            res += (permsCnt * partsCnt);
            permsCnt *= (i + 1);
            ResetP2D(p2d);
        }
    }
}

void CountCompsDistinctLen(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    mpz_class partsCnt = 0;
    mpz_class permsCnt = 0;

    NumPermsNoRepGmp(permsCnt, m, m);
    CountPartsDistinctLen(partsCnt, p1, p2, n, m);

    res = (partsCnt * permsCnt);
}


void CountCompsDistinctMultiZero(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        res = 1;
    } else {
        res = 0;
        mpz_class partsCnt = 0;
        mpz_class permsCnt = 0;
        NumPermsNoRepGmp(permsCnt, strtLen, strtLen);

        for (int i = strtLen; i <= m; ++i) {
            CountPartsDistinctLen(partsCnt, p1, p2, n, i);
            res += (permsCnt * partsCnt);
            permsCnt *= (i + 1);
        }
    }
}

void CountCompsDistinctMZWeak(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (strtLen == 0) {
        // This means that z contains only zeros
        res = 1;
    } else {
        mpz_class partsCnt = 1;
        mpz_class permsCnt = 1;

        for (int i = m; i > m - strtLen; --i) {
            permsCnt *= i;
        }

        res = 0;

        for (int i = strtLen; i <= m; ++i) {
            CountPartsDistinctLen(partsCnt, p1, p2, n, i);
            res += (permsCnt * partsCnt);
            permsCnt *= (m - i);
        }
    }
}
