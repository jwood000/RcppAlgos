#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountSection.h"
#include "Combinations/BigComboCount.h"
#include "Cpp14MakeUnique.h"
#include "SetUpUtils.h"

void CountPartsRepLenCap(mpz_t res, mpz_t* p1, mpz_t* p2,
                         int n, int m, int cap, int strtLen) {

    if (cap > n) cap = n;

    if (cap * m < n || n < m) {
        mpz_set_ui(res, 0u);
    } else if (cap * m == n || n <= m + 1) {
        mpz_set_ui(res, 1u);
    } else if (m < 2) {
        mpz_set_si(res, m);
    } else if (m == 2) {
        if (cap * 2 >= n) {
            cap = std::min(cap, n - 1);
            int result = n / m - (n - 1 - cap);
            mpz_set_ui(res, result);
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
                        mpz_add(p1[j + k], p2[j + k - 1], p1[j1 + k - i]);
                    }
                }
            } else {
                for (int i = 0; i < maxSize; ++i) {
                    mpz_set_ui(p2[i], 0u);
                }

                for (int j = width; j < maxSize; j += width) {
                    for (int k = i, j1 = j - width; k < width; ++k) {
                        mpz_add(p2[j + k], p1[j + k - 1], p2[j1 + k - i]);
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

void CountPartsRepLen(mpz_t res, mpz_t* p1, mpz_t* p2,
                      int n, int m, int cap, int strtLen) {

    if (m == 0 && n == 0) {
        mpz_set_ui(res, 1u);
    } else if (m == 0) {
        mpz_set_ui(res, 0u);
    } else if (n < m) {
        mpz_set_ui(res, 0u);
    } else if (n == m) {
        mpz_set_ui(res, 1u);
    } else if (m < 2) {
        mpz_set_ui(res, 1u);
    } else if (n - m == 1) {
        mpz_set_ui(res, 1u);
    } else if (m == 2) {
        int ndiv2 = n / 2;
        mpz_set_si(res, ndiv2);
    } else if (n - m == 2) {
        mpz_set_ui(res, 2u);
    } else if (m == 3) {
        mpz_t mpzN;
        mpz_init(mpzN);
        mpz_set_si(mpzN, n);
        SumSection(mpzN, res);
        mpz_clear(mpzN);
    } else {
        const int limit = std::min(n - m, m);
        n = (n < 2 * m) ? 2 * limit : n;

        if (n <= typeSwitchBnd) {
            for (int i = 3; i <= n; ++i) {
                double temp = SumSection(i);
                mpz_set_d(p1[i], temp);
            }
        } else {
            for (int i = 3; i < typeSwitchBnd; ++i) {
                double temp = SumSection(i);
                mpz_set_d(p1[i], temp);
            }

            mpz_t tempN;
            mpz_init(tempN);

            for (int i = typeSwitchBnd; i <= n; ++i) {
                mpz_set_si(tempN, i);
                SumSection(tempN, p1[i]);
            }

            mpz_clear(tempN);
        }

        for (int i = 4; i <= limit; ++i) {
            const int m2 = i * 2;

            if (i % 2) {
                mpz_set_ui(p1[i], 1u);

                for (int j = i + 1; j < m2; ++j) {
                    mpz_set(p1[j], p2[j - 1]);
                }

                for (int j = m2; j <= n; ++j) {
                    mpz_add(p1[j], p2[j - 1], p1[j - i]);
                }
            } else {
                mpz_set_ui(p2[i], 1u);

                for (int j = i + 1; j < m2; ++j) {
                    mpz_set(p2[j], p1[j - 1]);
                }

                for (int j = m2; j <= n; ++j) {
                    mpz_add(p2[j], p1[j - 1], p2[j - i]);
                }
            }
        }

        if (limit % 2) {
            mpz_set(res, p1[n]);
        } else {
            mpz_set(res, p2[n]);
        }
    }
}

void CountPartsRep(mpz_t res, int n, int m, int cap, int strtLen) {

    auto qq = FromCpp14::make_unique<mpz_t[]>(n + 1);

    for (int i = 0; i <= n; ++i) {
        mpz_init(qq[i]);
    }

    mpz_set_ui(qq[0], 1u);
    mpz_set_ui(qq[1], 1u);

    for(int i = 2; i <= n; ++i) {
        mpz_set_ui(qq[i], 0u);

        for (int s = 1, f = 1, r = 1; i >= r; f += 3, r += f, s *= -1) {
            if (s > 0) {
                mpz_add(qq[i], qq[i], qq[i - r]);
            } else {
                mpz_sub(qq[i], qq[i], qq[i - r]);
            }
        }

        for (int s = 1, f = 2, r = 2; i >= r; f += 3, r += f, s *= -1) {
            if (s > 0) {
                mpz_add(qq[i], qq[i], qq[i - r]);
            } else {
                mpz_sub(qq[i], qq[i], qq[i - r]);
            }
        }
    }

    mpz_set(res, qq[n]);
    MpzClearVec(qq.get(), n + 1);
}

void CountCompsRepLen(mpz_t res, int n, int m, int cap, int strtLen) {
    nChooseKGmp(res, n - 1, m - 1);
}

void CountCompsRepZero(mpz_t res, int n, int m, int cap, int strtLen) {

    if (n == m) {
        mpz_set_ui(res, 1);
        mpz_mul_2exp(res, res, n - 1);
    } else {
        mpz_t temp;
        mpz_init(temp);

        for (int i = 1; i <= m; ++i) {
            nChooseKGmp(temp, n - 1, i - 1);
            mpz_add(res, res, temp);
        }

        mpz_clear(temp);
    }
}
