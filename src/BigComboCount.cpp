#include "Cpp14MakeUnique.h"
#include <vector>
#include <algorithm> // std::min
#include <numeric>   // std::accumulate
#include <deque>
#include <gmp.h>

// All functions below are exactly the same as the functions
// in StdCombinationCount.cpp. The only difference is that they
// utilize the gmp library and deal mostly with mpz_t types

void nChooseKGmp(mpz_t result, int n, int k) {
    mpz_bin_uiui(result, n, k);
}

void NumCombsWithRepGmp(mpz_t result, int n, int k) {
    mpz_bin_uiui(result, n + k - 1, k);
}

bool OnlyOneCombo(int n, int k, const std::deque<int> &Reps) {
    if (k < 1 || n <= 1) {
        return true;
    }

    if (k == n && std::accumulate(Reps.cbegin(),
                                  Reps.cend(), 0) == n) {
        return true;
    }

    return false;
}

void MultisetCombRowNumGmp(mpz_t result, int n, int r,
                           const std::deque<int> &Reps) {

    if (!OnlyOneCombo(n, r, Reps)) {
        const int r1 = r + 1;
        int myMax = r1;

        if (myMax > Reps[0] + 1) {
            myMax = Reps[0] + 1;
        }

        auto triangleVec = FromCpp14::make_unique<mpz_t[]>(r1);
        auto temp = FromCpp14::make_unique<mpz_t[]>(r1);

        for (int i = 0; i < r1; ++i) {
            mpz_init(triangleVec[i]);
            mpz_init(temp[i]);
        }

        for (int i = 0; i < myMax; ++i) {
            mpz_set_ui(triangleVec[i], 1);
            mpz_set_ui(temp[i], 1);
        }

        --myMax;
        int ind = 1;

        for (; myMax < r; ++ind) {
            int myMin = std::min(Reps[ind], r);

            for (int i = 1; i <= myMin; ++i) {
                mpz_add(triangleVec[i], triangleVec[i], triangleVec[i - 1]);
            }

            myMin = std::min(Reps[ind] + myMax, r);
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
            mpz_set(t, triangleVec[r]);
            const int s = std::min(Reps[ind], r);

            for (int i = 1; i <= s; ++i) {
                mpz_add(triangleVec[r], triangleVec[r],
                        triangleVec[r - i]);
            }

            mpz_set(mySum, triangleVec[r]);

            for (int i = r - 1; i >= s; --i) {
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
            const int myMin2 = std::min(Reps[n1], r);

            for (int i = 1; i <= myMin2; ++i) {
                mpz_add(triangleVec[r], triangleVec[r],
                        triangleVec[r - i]);
            }
        }

        mpz_set(result, triangleVec[r]);

        for (int i = 0; i < r1; ++i) {
            mpz_clear(triangleVec[i]);
            mpz_clear(temp[i]);
        }

        mpz_clear(mySum);
        mpz_clear(t);
    } else {
        mpz_set_ui(result, 1);
    }
}
