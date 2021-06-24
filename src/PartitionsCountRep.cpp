#include "Partitions/PartitionsCountSection.h"
#include "Combinations/ComboCount.h"
#include "Cpp14MakeUnique.h"
#include "CleanConvert.h"
#include <vector>
#include <cmath>

int width;
int blockSize;
static std::vector<double> memoize;

double pStdCap(int n, int m, int cap) {

    if (cap * m < n || n < m) return 0;
    if (cap * m == n || n <= m + 1) return 1;
    if (m < 2) return m;

    const int block = cap * blockSize + (n - m) * width + m - 2;
    if (memoize[block]) return memoize[block];

    int niter = n / m;

    if (m == 2) {
        if (cap * 2 >= n) {
            cap = std::min(cap, n - 1);
            return niter - (n - 1 - cap);
        } else {
            return 0;
        }
    }

    double count = 0;

    for (; niter--; n -= m, --cap) {
        count += (memoize[cap * blockSize + (n - m) * width + m - 3] =
                  pStdCap(n - 1, m - 1, cap));
    }

    return count;
}

double CountPartRepLenCap(int n, int m, int cap) {

    if (cap > n) cap = n;
    if (cap * m < n || n < m) return 0;
    if (cap * m == n || n <= m + 1) return 1;
    if (m < 2) return m;

    if (m == 2) {
        if (cap * 2 >= n) {
            cap = std::min(cap, n - 1);
            return n / m - (n - 1 - cap);
        } else {
            return 0;
        }
    }

    width = m;
    blockSize = m * (n - m + 1);
    memoize = std::vector<double>((cap + 1) * blockSize, 0.0);
    return pStdCap(n, m, cap);
}

double CountPartRepLen(int n, int m) {

    if (m == 0)
        return (n == 0) ? 1.0 : 0.0;

    if (n < m)
        return 0.0;

    if (n == m)
        return 1.0;

    if (m < 2)
        return 1.0;

    if (n - m == 1)
        return 1.0;

    // If n > 3, we have the following:
    // 1st part: n - 3 1's followed by a 3
    // 2nd part: n - 4 1's followed by two 2's
    //
    // N.B. We have taken care of every case where n <= 2 above
    // i.e. n = 2, m = 2; n = 2, m = 1; n = 1, m = 1; n = 0, m = 0
    if (n - m == 2)
        return 2.0;

    if (m == 2)
        return std::floor(n / 2);

    if (m == 3) {
        const double res = SumSection(n);
        return(res);
    } else {
        const int limit = std::min(n - m, m);
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

void CountPartRepLen(mpz_t res, int n, int m) {

    if (m == 3) {
        mpz_t mpzN;
        mpz_init(mpzN);
        mpz_set_si(mpzN, n);
        SumSection(mpzN, res);
        mpz_clear(mpzN);
    } else {
        const int limit = std::min(n - m, m);
        n = (n < 2 * m) ? 2 * limit : n;

        auto p1 = FromCpp14::make_unique<mpz_t[]>(n + 1);
        auto p2 = FromCpp14::make_unique<mpz_t[]>(n + 1);

        for (int i = 0; i <= n; ++i) {
            mpz_init(p1[i]);
            mpz_init(p2[i]);
        }

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

        for (int i = 0; i <= n; ++i) {
            mpz_clear(p1[i]);
            mpz_clear(p2[i]);
        }
    }
}

// Similar to CountPartDistinct
double CountPartRep(int n) {

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

void CountPartRep(mpz_t res, int n) {

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

    for (int i = 0; i <= n; ++i) {
        mpz_clear(qq[i]);
    }
}

double CountPartPermRep(int target, int m, bool includeZero) {
    return (includeZero) ? nChooseK(target + m - 1, m - 1) :
                           nChooseK(target - 1, m - 1);
}
