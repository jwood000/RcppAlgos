#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountSection.h"
#include "Cpp14MakeUnique.h"
#include <vector>

void CountPartsDistinctLenCap(mpz_t res, mpz_t* p1, mpz_t* p2,
                              int n, int m, int cap, int strtLen) {

    if (cap > n) cap = n;
    const int limit = (cap * m) - ((m - 1) * m) / 2;

    if (m > n || cap < m) {
        mpz_set_ui(res, 0u);
    } else if (m == n) {
        if (n == 1 && cap >= 1) {
            mpz_set_ui(res, 1u);
        } else {
            mpz_set_ui(res, 0u);
        }
    } else if (m == 1) {
        if (cap >= n) {
            mpz_set_ui(res, 1u);
        } else {
            mpz_set_ui(res, 0u);
        }
    } else if (limit <= n) {
        if (limit == n) {
            mpz_set_ui(res, 1u);
        } else {
            mpz_set_ui(res, 0u);
        }
    } else {
        const int width = n + 1;
        const int maxSize = (cap + 1) * width;

        for (int i = 0; i < maxSize; ++i) {
            mpz_set_ui(p1[i], 0u);
        }

        for (int i = 1; i < width; ++i) {
            for (int j = i; j <= cap; ++j) {
                mpz_set_ui(p1[j * width + i], 1u);
            }
        }

        for (int i = 2; i <= m; ++i) {
            if (i % 2) {
                for (int i = 0; i < maxSize; ++i) {
                    mpz_set_ui(p1[i], 0u);
                }

                for (int j = width; j < maxSize; j += width) {
                    for (int k = i, j1 = j - width; k < width; ++k) {
                        mpz_add(p1[j + k], p2[j1 + k - i], p1[j1 + k - i]);
                    }
                }
            } else {
                for (int i = 0; i < maxSize; ++i) {
                    mpz_set_ui(p2[i], 0u);
                }

                for (int j = width; j < maxSize; j += width) {
                    for (int k = i, j1 = j - width; k < width; ++k) {
                        mpz_add(p2[j + k], p1[j1 + k - i], p2[j1 + k - i]);
                    }
                }
            }
        }

        if (m % 2) {
            mpz_set(res, p1[maxSize - 1]);
        } else {
            mpz_set(res, p2[maxSize - 1]);
        }
    }
}

void CountPartsDistinctLen(mpz_t res, mpz_t* p1, mpz_t* p2,
                           int n, int m, int cap, int strtLen) {

    const int max_width = GetMaxWidth(n);

    if (m == 0 && n == 0) {
        mpz_set_ui(res, 1u);
    } else if (m == 0) {
        mpz_set_ui(res, 0u);
    } else if (m > max_width) {
        mpz_set_ui(res, 0u);
    } else if (m < 2) {
        mpz_set_ui(res, 1u);
    } else if (m == 2) {
        int n1div2 = (n - 1) / 2;
        mpz_set_si(res, n1div2);
    } else if (m == 3) {
        mpz_t mpzN;
        mpz_init(mpzN);
        mpz_set_si(mpzN, n - 3);
        SumSection(mpzN, res);
        mpz_clear(mpzN);
    } else {
        const int limit = (m == GetMaxWidth(n + 1)) ? m - 1 : m;

        if (n <= typeSwitchBnd) {
            for (int i = 6; i <= n; ++i) {
                double temp = SumSection(i - 3);
                mpz_set_d(p1[i], temp);
            }
        } else {
            for (int i = 6; i < typeSwitchBnd; ++i) {
                double temp = SumSection(i - 3);
                mpz_set_d(p1[i], temp);
            }

            mpz_t tempN;
            mpz_init(tempN);

            for (int i = typeSwitchBnd; i <= n; ++i) {
                mpz_set_si(tempN, i - 3);
                SumSection(tempN, p1[i]);
            }
        }

        for (int i = 4; i <= limit; ++i) {
            const int m1 = ((i + 1) * i) / 2;
            const int m2 = m1 + i;

            if (i % 2) {
                for (int j = m1; j < m2; ++j) {
                    mpz_set(p1[j], p2[j - i]);
                }

                for (int j = m2; j <= n; ++j) {
                    mpz_add(p1[j], p2[j - i], p1[j - i]);
                }
            } else {
                for (int j = m1; j < m2; ++j) {
                    mpz_set(p2[j], p1[j - i]);
                }

                for (int j = m2; j <= n; ++j) {
                    mpz_add(p2[j], p2[j - i], p1[j - i]);
                }
            }
        }

        if (m > limit && m % 2) {
            mpz_set(res, p2[n - m]);
        } else if (m > limit) {
            mpz_set(res, p1[n - m]);
        } else if (m % 2) {
            mpz_set(res, p1[n]);
        } else {
            mpz_set(res, p2[n]);
        }
    }
}

void CountPartsDistinct(mpz_t res, int n, int m, int cap, int strtLen) {

    auto qq = FromCpp14::make_unique<mpz_t[]>(n + 1);

    mpz_init(qq[0]);
    mpz_init(qq[1]);

    mpz_set_ui(qq[0], 1u);
    mpz_set_ui(qq[1], 1u);

    for(int i = 2 ; i <= n; ++i) {
        mpz_init(qq[i]);

        for (int s = 1, f = 5, r = 2; i >= r; r += f, f += 3, s *= -1) {
            if (s > 0) {
                mpz_add(qq[i], qq[i], qq[i - r]);

                if (i == r * 2) {
                    mpz_sub_ui(qq[i], qq[i], 1u);
                }
            } else {
                mpz_sub(qq[i], qq[i], qq[i - r]);

                if (i == r * 2) {
                    mpz_add_ui(qq[i], qq[i], 1u);
                }
            }
        }

        for (int s = 1, f = 4, r = 1; i >= r; r += f, f += 3, s *= -1) {
            if (s > 0) {
                mpz_add(qq[i], qq[i], qq[i - r]);

                if (i == r * 2) {
                    mpz_sub_ui(qq[i], qq[i], 1u);
                }
            } else {
                mpz_sub(qq[i], qq[i], qq[i - r]);

                if (i == r * 2) {
                    mpz_add_ui(qq[i], qq[i], 1u);
                }
            }
        }
    }

    mpz_set(res, qq[n]);

    for (int i = 0; i <= n; ++i) {
        mpz_clear(qq[i]);
    }
}

void CountPartsDistinctMultiZero(mpz_t res, mpz_t* p1, mpz_t* p2,
                                 int n, int m, int cap, int strtLen) {

    mpz_t temp;
    mpz_init(temp);
    mpz_set_ui(res, 0);

    for (int i = strtLen; i <= m; ++i) {
        CountPartsDistinctLen(temp, p1, p2, n, i, cap, strtLen);
        mpz_add(res, res, temp);
    }

    mpz_clear(temp);
}

void CountPartsDistinctCapMZ(mpz_t res, mpz_t* p1, mpz_t* p2,
                             int n, int m, int cap, int strtLen) {

    mpz_t temp;
    mpz_init(temp);
    mpz_set_ui(res, 0);

    for (int i = strtLen; i <= m; ++i) {
        CountPartsDistinctLenCap(temp, p1, p2, n, i, cap, strtLen);
        mpz_add(res, res, temp);
    }

    mpz_clear(temp);
}
