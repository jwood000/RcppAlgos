#include <algorithm> // std::min, std::fill_n
#include "CppConvert/GmpxxCopy.h"
#include <numeric>   // std::accumulate
#include <vector>
#include <deque>

// All functions below are exactly the same as the functions
// in StdCombinationCount.cpp. The only difference is that they
// utilize the gmp library and deal mostly with mpz_class types

void nChooseKGmp(mpz_class &result, int n, int m) {
    mpz_bin_uiui(result.get_mpz_t(), n, m);
}

void NumCombsWithRepGmp(mpz_class &result, int n, int m) {
    mpz_bin_uiui(result.get_mpz_t(), n + m - 1, m);
}

bool OnlyOneCombo(int n, int m, const std::deque<int> &Reps) {
    if (m < 1 || n <= 1) return true;
    return (m == n && std::accumulate(Reps.cbegin(), Reps.cend(), 0) == n);
}

void MultisetCombRowNumGmp(mpz_class &result, int n, int m,
                           const std::deque<int> &Reps) {

    if (!OnlyOneCombo(n, m, Reps)) {
        const int m1 = m + 1;
        int myMax = m1;

        if (myMax > Reps[0] + 1) {
            myMax = Reps[0] + 1;
        }

        std::vector<mpz_class> triangleVec(m1);
        std::vector<mpz_class> temp(m1);

        std::fill_n(triangleVec.begin(), myMax, 1);
        std::fill_n(temp.begin(), myMax, 1);

        --myMax;
        int ind = 1;

        for (; myMax < m; ++ind) {
            int myMin = std::min(Reps[ind], m);

            for (int i = 1; i <= myMin; ++i) {
                triangleVec[i] += triangleVec[i - 1];
            }

            myMin = std::min(Reps[ind] + myMax, m);
            int j = 0;

            for (int i = (Reps[ind] + 1); i <= myMin; ++i, ++j) {
                triangleVec[i] += triangleVec[i - 1];
                triangleVec[i] -= temp[j];
                temp[j] = triangleVec[j];
            }

            for (; j <= myMin; ++j) {
                temp[j] = triangleVec[j];
            }

            myMax = myMin;
        }

        const int n1 = n - 1;

        mpz_class t;
        mpz_class mySum;

        for (; ind < n1; ++ind) {
            t = triangleVec[m];
            const int s = std::min(Reps[ind], m);

            for (int i = 1; i <= s; ++i) {
                triangleVec[m] += triangleVec[m - i];
            }

            mySum = triangleVec[m];

            for (int i = m - 1; i >= s; --i) {
                mySum -= t;
                t = triangleVec[i];
                mySum += triangleVec[i - s];
                triangleVec[i] = mySum;
            }

            for (int i = s - 1; i > 0; --i) {
                mySum -= t;
                t = triangleVec[i];
                triangleVec[i] = mySum;
            }
        }

        if (ind < n) {
            const int myMin2 = std::min(Reps[n1], m);

            for (int i = 1; i <= myMin2; ++i) {
                triangleVec[m] += triangleVec[m - i];
            }
        }

        result = triangleVec[m];
    } else {
        result = 1;
    }
}
