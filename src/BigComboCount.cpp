#include "Cpp14MakeUnique.h"
#include <vector>
#include <algorithm> // std::min
#include <numeric>   // std::accumulate
#include <deque>
#include <gmp.h>

// All functions below are exactly the same as the functions
// in StdCombinationCount.cpp. The only difference is that they
// utilize the gmp library and deal mostly with mpz_t types

void nChooseKGmp(mpz_t result, int n, int m) {
    mpz_bin_uiui(result, n, m);
}

void NumCombsWithRepGmp(mpz_t result, int n, int m) {
    mpz_bin_uiui(result, n + m - 1, m);
}

bool OnlyOneCombo(int n, int m, const std::deque<int> &Reps) {
    if (m < 1 || n <= 1) {
        return true;
    }

    if (m == n && std::accumulate(Reps.cbegin(),
                                  Reps.cend(), 0) == n) {
        return true;
    }

    return false;
}

void MultisetCombRowNumGmp(mpz_t result, int n, int m,
                           const std::deque<int> &Reps) {

    if (!OnlyOneCombo(n, m, Reps)) {
        const int m1 = m + 1;
        int myMax = m1;

        if (myMax > Reps[0] + 1) {
            myMax = Reps[0] + 1;
        }

        auto triangleVec = FromCpp14::make_unique<mpz_t[]>(m1);
        auto temp = FromCpp14::make_unique<mpz_t[]>(m1);

        for (int i = 0; i < m1; ++i) {
            mpz_init(triangleVec[i]);
            mpz_init(temp[i]);
        }

        for (int i = 0; i < myMax; ++i) {
            mpz_set_ui(triangleVec[i], 1);
            mpz_set_ui(temp[i], 1);
        }

        --myMax;
        int ind = 1;

        for (; myMax < m; ++ind) {
            int myMin = std::min(Reps[ind], m);

            for (int i = 1; i <= myMin; ++i) {
                mpz_add(triangleVec[i], triangleVec[i], triangleVec[i - 1]);
            }

            myMin = std::min(Reps[ind] + myMax, m);
            int j = 0;

            for (int i = (Reps[ind] + 1); i <= myMin; ++i, ++j) {
                mpz_add(triangleVec[i], triangleVec[i], triangleVec[i - 1]);
                mpz_sub(triangleVec[i], triangleVec[i], temp[j]);
                mpz_set(temp[j], triangleVec[j]);
            }

            for (; j <= myMin; ++j) {
                mpz_set(temp[j], triangleVec[j]);
            }

            myMax = myMin;
        }

        const int n1 = n - 1;

        mpz_t t;
        mpz_t mySum;

        mpz_init(t);
        mpz_init(mySum);

        for (; ind < n1; ++ind) {
            mpz_set(t, triangleVec[m]);
            const int s = std::min(Reps[ind], m);

            for (int i = 1; i <= s; ++i) {
                mpz_add(triangleVec[m], triangleVec[m],
                        triangleVec[m - i]);
            }

            mpz_set(mySum, triangleVec[m]);

            for (int i = m - 1; i >= s; --i) {
                mpz_sub(mySum, mySum, t);
                mpz_set(t, triangleVec[i]);
                mpz_add(mySum, mySum, triangleVec[i - s]);
                mpz_set(triangleVec[i], mySum);
            }

            for (int i = s - 1; i > 0; --i) {
                mpz_sub(mySum, mySum, t);
                mpz_set(t, triangleVec[i]);
                mpz_set(triangleVec[i], mySum);
            }
        }

        if (ind < n) {
            const int myMin2 = std::min(Reps[n1], m);

            for (int i = 1; i <= myMin2; ++i) {
                mpz_add(triangleVec[m], triangleVec[m],
                        triangleVec[m - i]);
            }
        }

        mpz_set(result, triangleVec[m]);

        for (int i = 0; i < m1; ++i) {
            mpz_clear(triangleVec[i]);
            mpz_clear(temp[i]);
        }

        mpz_clear(mySum);
        mpz_clear(t);
    } else {
        mpz_set_ui(result, 1);
    }
}
