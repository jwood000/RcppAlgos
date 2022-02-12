#include "Permutations/PermuteCount.h"
#include "Cpp14MakeUnique.h"
#include <algorithm> // std::sort, std::max_element
#include <numeric>   // std::accumulate, std::iota
#include <gmp.h>

// All functions below are exactly the same as the functions
// in StandardCount.cpp. The only difference is that they
// utilize the gmp library and deal mostly with mpz_t types

void NumPermsWithRepGmp(mpz_t result, const std::vector<int> &v) {

    mpz_set_ui(result, 1u);
    std::vector<int> myLens = rleCpp(v);
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());

    const int myMax = myLens[0];
    const int numUni = myLens.size();

    for (int i = v.size(); i > myMax; --i)
        mpz_mul_ui(result, result, i);

    if (numUni > 1) {
        mpz_t div;
        mpz_init_set_ui(div, 1u);

        for (int i = 1; i < numUni; ++i) {
            for (int j = 2; j <= myLens[i]; ++j) {
                mpz_mul_ui(div, div, j);
            }
        }

        mpz_divexact(result, result, div);
        mpz_clear(div);
    }
}

void NumPermsNoRepGmp(mpz_t result, int n, int k) {

    mpz_set_ui(result, 1u);

    for (int i = n, m = n - k; i > m; --i)
        mpz_mul_ui(result, result, i);
}

void MultisetPermRowNumGmp(mpz_t result, int n, int r,
                           const std::vector<int> &myReps) {

    const int sumFreqs = std::accumulate(myReps.cbegin(), myReps.cend(), 0);

    if (n < 2 || r < 1) {
        mpz_set_ui(result, 1);
    } else if (r > sumFreqs) {
        mpz_set_ui(result, 0);
    } else if (r == sumFreqs) {
        std::vector<int> freqs(sumFreqs);

        for (int i = 0, k = 0; i < static_cast<int>(myReps.size()); ++i)
            for (int j = 0; j < myReps[i]; ++j, ++k)
                freqs[k] = i;

        NumPermsWithRepGmp(result, freqs);
    } else {
        const int n1 = n - 1;
        int maxFreq = *std::max_element(myReps.cbegin(), myReps.cend());

        std::vector<int> seqR(r);
        std::iota(seqR.begin(), seqR.end(), 1);

        mpz_t prodR;
        mpz_init(prodR);
        mpz_set_ui(prodR, 1);

        for (int i = 0; i < r; ++i) {
            mpz_mul_ui(prodR, prodR, seqR[i]);
        }

        const std::size_t uR1 = r + 1;
        const int myMax = (r < maxFreq) ? (r + 2) : (maxFreq + 2);
        auto cumProd = FromCpp14::make_unique<mpz_t[]>(myMax);
        auto resV = FromCpp14::make_unique<mpz_t[]>(uR1);

        for (int i = 0; i < myMax; ++i) {
            mpz_init(cumProd[i]);
        }

        // Equivalent to c(1, 1:myMax)
        mpz_set_ui(cumProd[0], 1);

        for (int i = 1; i < myMax; ++i) {
            mpz_set_ui(cumProd[i], i);
        }

        for (std::size_t i = 0; i < uR1; ++i) {
            mpz_init(resV[i]);
            mpz_set_ui(resV[i], 0);
        }

        for (int i = 1; i < myMax; ++i) {
            mpz_mul(cumProd[i], cumProd[i], cumProd[i - 1]);
        }

        int myMin = std::min(r, myReps[0]);

        for (int i = 0; i <= myMin; ++i) {
            mpz_divexact(resV[i], prodR, cumProd[i]);
        }

        mpz_t temp;
        mpz_init(temp);

        for (int i = 1; i < n1; ++i) {
            for (int j = r; j > 0; --j) {
                myMin = std::min(j, myReps[i]);
                mpz_set_ui(result, 0);

                for (int k = 0; k <= myMin; ++k) {
                    mpz_divexact(temp, resV[j - k], cumProd[k]);
                    mpz_add(result, result, temp);
                }

                mpz_set(resV[j], result);
            }
        }

        myMin = std::min(r, myReps[n1]);
        mpz_set_ui(result, 0);

        for (int k = 0; k <= myMin; ++k) {
            mpz_divexact(temp, resV[r - k], cumProd[k]);
            mpz_add(result, result, temp);
        }

        mpz_clear(temp);
        mpz_clear(prodR);

        for (int i = 0; i < myMax; ++i) {
            mpz_clear(cumProd[i]);
        }

        for (std::size_t i = 0; i < uR1; ++i) {
            mpz_clear(resV[i]);
        }
    }
}
