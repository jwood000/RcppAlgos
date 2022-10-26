#include "Permutations/PermuteCount.h"
#include <algorithm> // std::sort, std::max_element
#include <numeric>   // std::accumulate, std::iota
#include <gmpxx.h>

// All functions below are exactly the same as the functions
// in StandardCount.cpp. The only difference is that they
// utilize the gmp library and deal mostly with mpz_t types

void NumPermsWithRepGmp(mpz_class &result, const std::vector<int> &v) {

    result = 1;
    std::vector<int> myLens = rleCpp(v);
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());

    const int myMax = myLens[0];
    const int numUni = myLens.size();

    for (int i = v.size(); i > myMax; --i) {
        result *= i;
    }

    if (numUni > 1) {
        mpz_class div(1);

        for (int i = 1; i < numUni; ++i) {
            for (int j = 2; j <= myLens[i]; ++j) {
                div *= j;
            }
        }

        mpz_divexact(result.get_mpz_t(), result.get_mpz_t(), div.get_mpz_t());
    }
}

void NumPermsNoRepGmp(mpz_class &result, int n, int k) {

    result = 1;

    for (int i = n, m = n - k; i > m; --i) {
        result *= i;
    }
}

void MultisetPermRowNumGmp(mpz_class &result, int n, int m,
                           const std::vector<int> &myReps) {

    const int sumFreqs = std::accumulate(myReps.cbegin(), myReps.cend(), 0);

    if (n < 2 || m < 1) {
        result = 1;
    } else if (m > sumFreqs) {
        result = 0;
    } else if (m == sumFreqs) {
        std::vector<int> freqs(sumFreqs);

        for (int i = 0, k = 0; i < static_cast<int>(myReps.size()); ++i) {
            for (int j = 0; j < myReps[i]; ++j, ++k) {
                freqs[k] = i;
            }
        }

        NumPermsWithRepGmp(result, freqs);
    } else {
        const int n1 = n - 1;
        int maxFreq = *std::max_element(myReps.cbegin(), myReps.cend());

        std::vector<int> seqR(m);
        std::iota(seqR.begin(), seqR.end(), 1);

        mpz_class prodR(1);

        for (int i = 0; i < m; ++i) {
            prodR *= seqR[i];
        }

        const std::size_t uR1 = m + 1;
        const int myMax = (m < maxFreq) ? (m + 2) : (maxFreq + 2);

        std::vector<mpz_class> cumProd(myMax);
        std::vector<mpz_class> resV(uR1, 0);

        // Equivalent to c(1, 1:myMax)
        cumProd[0] = 1;

        for (int i = 1; i < myMax; ++i) {
            cumProd[i] = i;
        }

        for (int i = 1; i < myMax; ++i) {
            cumProd[i] *= cumProd[i - 1];
        }

        int myMin = std::min(m, myReps[0]);

        for (int i = 0; i <= myMin; ++i) {
            mpz_divexact(resV[i].get_mpz_t(), prodR.get_mpz_t(),
                         cumProd[i].get_mpz_t());
        }

        mpz_class temp;

        for (int i = 1; i < n1; ++i) {
            for (int j = m; j > 0; --j) {
                myMin = std::min(j, myReps[i]);
                result = 0;

                for (int k = 0; k <= myMin; ++k) {
                    mpz_divexact(temp.get_mpz_t(), resV[j - k].get_mpz_t(),
                                 cumProd[k].get_mpz_t());
                    result += temp;
                }

                resV[j] = result;
            }
        }

        myMin = std::min(m, myReps[n1]);
        result = 0;

        for (int k = 0; k <= myMin; ++k) {
            mpz_divexact(temp.get_mpz_t(), resV[m - k].get_mpz_t(),
                         cumProd[k].get_mpz_t());
            result += temp;
        }
    }
}
