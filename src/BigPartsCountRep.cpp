#include "Partitions/PartitionsCountSection.h"
#include "Partitions/BigPartsCountSection.h"
#include "Combinations/BigComboCount.h"
#include <algorithm>
#include <vector>

// See commentary in PartitionsCountRep.cpp
void CountPartsRepLenRstrctd(
    mpz_class &res, std::vector<std::vector<mpz_class>> &p2d,
    int n, int m, const std::vector<int> &allowed, int strtLen
) {

    CheckAllowedInvariant(allowed);
    ResetP2D(p2d);
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

// *********************** CountCompsRepLenCap Details ***********************
//
// This uses the inclusion–exclusion principle to count the number of positive
// compositions of length m that sum to n, with each part bounded above by cap.
//
// First, count all positive compositions without the cap:
//
//     x_1 + x_2 + ... + x_m = n,    x_i >= 1
//
// This is:
//
//     choose(n - 1, m - 1)
//
// Next, subtract compositions where at least one part exceeds the cap.
// Suppose exactly one part violates the cap. Without loss of generality:
//
//     x_1 >= cap + 1
//
// Write:
//
//     x_1 = cap + z_1
//     x_i = z_i   for i >= 2
//
// Then:
//
//     z_1 + z_2 + ... + z_m = n - cap
//
// The number of such compositions is:
//
//     choose(n - cap - 1, m - 1)
//
// and there are choose(m, 1) choices for which part violates the cap.
//
// However, this subtracts too much, because compositions where two parts
// exceed the cap were subtracted twice. To correct this, add back the number
// of compositions where two parts exceed the cap:
//
//     z_1 + ... + z_m = n - 2*cap
//
// giving:
//
//     choose(n - 2*cap - 1, m - 1)
//
// and choose(m, 2) ways to select the violating parts.
//
// Continuing in this way, alternating subtraction and addition, yields:
//
//     sum_{i=0}^{floor((n - m)/cap)}
//         (-1)^i * choose(m, i) * choose(n - i*cap - 1, m - 1)
//
// The upper limit floor((n - m)/cap) arises because each violating part must
// contribute at least cap extra beyond the minimum value of 1.
//
// Count capped positive compositions with repetition.
// PRECONDITIONS (enforced by caller / dispatch):
//   - n >= m >= 1
//   - cap >= 2
//   - n <= cap * m
//   - allowed is a single-element vector containing cap (allowed.size() == 1),
//     or the set [1..cap] defined in PartitionsCount
//
// NOTE: This function does NOT support arbitrary allowed-sets. It counts
//       parts in the standard range [1..cap].
void CountCompsRepLenCap(mpz_class &res, int n, int m,
                         const std::vector<int> &allowed, int strtLen = 0) {

    static_cast<void>(strtLen); // intentionally unused

    const int cap = *std::max_element(allowed.cbegin(), allowed.cend());
    const int maxViolations = (n - m) >= 0 ? (n - m) / cap : -1;

    mpz_class nSlots;
    mpz_class badParts;
    nChooseKGmp(res, n - 1, m - 1);

    for (int i = 1, sign = -1; i <= maxViolations; ++i, sign *= -1) {
        nChooseKGmp(nSlots, m, i);
        nChooseKGmp(badParts, n - i * cap - 1, m - 1);
        res += (sign * nSlots * badParts);
    }
}

void CountCompsRepCapZNotWk(mpz_class &res, int n, int m,
                            const std::vector<int> &allowed, int strtLen) {

    static_cast<void>(strtLen);

    res = 0;
    mpz_class temp;

    const int cap = *std::max_element(allowed.cbegin(), allowed.cend());
    const int minWidth = ((n - 1) / cap) + 1;  // ceil(n / cap)

    for (int i = minWidth; i <= m; ++i) {
        temp = 0; // defensive: avoid accumulation if the callee doesn't reset
        CountCompsRepLenCap(temp, n, i, allowed);
        res += temp;
    }
}

void CountCompsRepZNotWk(mpz_class &res, int n, int m,
                         const std::vector<int> &allowed, int strtLen) {

    static_cast<void>(allowed);
    static_cast<void>(strtLen);
    res = 0;

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
