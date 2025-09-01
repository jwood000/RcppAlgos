#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountSection.h"
#include "Combinations/BigComboCount.h"
#include <vector>

void CountPartsRepLenRstrctd(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    p2d[0][0] = 1; // one way to make 0 with 0 parts

    for (int num : allowed) {
        if (num) {
            for (int j = 1; j <= m; ++j) {
                for (int s = num; s <= n; ++s) {
                    p2d[j][s] += p2d[j - 1][s - num];
                }
            }
        }
    }

    res = p2d[m][n];
}

void CountPartsRepLen(
    mpz_class &res, std::vector<mpz_class> &p1, std::vector<mpz_class> &p2,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    if (m == 0 && n == 0) {
        res = 1;
    } else if (m == 0) {
        res = 0;
    } else if (n < m) {
        res = 0;
    } else if (n == m) {
        res = 1;
    } else if (m < 2) {
        res = 1;
    } else if (n - m == 1) {
        res = 1;
    } else if (m == 2) {
        res = (n / 2);
    } else if (n - m == 2) {
        res = 2;
    } else if (m == 3) {
        mpz_class mpzN(n);
        SumSection(mpzN, res);
    } else {
        const int limit = std::min(n - m, m);
        n = (n < 2 * m) ? 2 * limit : n;

        if (n <= typeSwitchBnd) {
            for (int i = 3; i <= n; ++i) {
                p1[i] = static_cast<double>(SumSection(i));
            }
        } else {
            for (int i = 3; i < typeSwitchBnd; ++i) {
                p1[i] = static_cast<double>(SumSection(i));
            }

            mpz_class tempN;

            for (int i = typeSwitchBnd; i <= n; ++i) {
                tempN = i;
                SumSection(tempN, p1[i]);
            }
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

        if (limit % 2) {
            res = p1[n];
        } else {
            res = p2[n];
        }
    }
}

void CountPartsRep(mpz_class &res, int n, int m,
                   const std::vector<int> &allowed, int strtLen) {

    std::vector<mpz_class> qq(n + 1, 0);

    qq[0] = 1;
    qq[1] = 1;

    for(int i = 2; i <= n; ++i) {
        for (int s = 1, f = 1, r = 1; i >= r; f += 3, r += f, s *= -1) {
            if (s > 0) {
                qq[i] += qq[i - r];
            } else {
                qq[i] -= qq[i - r];
            }
        }

        for (int s = 1, f = 2, r = 2; i >= r; f += 3, r += f, s *= -1) {
            if (s > 0) {
                qq[i] += qq[i - r];
            } else {
                qq[i] -= qq[i - r];
            }
        }
    }

    res = qq[n];
}

void CountCompsRepLen(mpz_class &res, int n, int m,
                      const std::vector<int> &allowed, int strtLen) {
    nChooseKGmp(res, n - 1, m - 1);
}

void CountCompsRepZNotWk(mpz_class &res, int n, int m,
                         const std::vector<int> &allowed, int strtLen) {

    if (n == m) {
        res = 1;
        mpz_mul_2exp(res.get_mpz_t(), res.get_mpz_t(), n - 1);
    } else {
        mpz_class temp;

        for (int i = 1; i <= m; ++i) {
            nChooseKGmp(temp, n - 1, i - 1);
            res += temp;
        }
    }
}
