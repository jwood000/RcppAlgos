#include "Permutations/BigPermuteCount.h"
#include "CleanConvert.h"
#include <algorithm>  // std::min, std::max, std::sort
#include <numeric>    // std::accumulate, std::partial_sum, std::iota
#include <limits>     // std::numeric_limits

// Most of the code for rleCpp was obtained from Hadley Wickham's
// article titled "High Performance functions with Rcpp" found:
//             http://adv-r.had.co.nz/Rcpp.html
std::vector<int> rleCpp(const std::vector<int> &x) {
    std::vector<int> lengths;
    int prev = x[0];
    std::size_t i = 0;
    lengths.push_back(1);

    for(auto it = x.cbegin() + 1; it != x.cend(); ++it) {
        if (prev == *it) {
            ++lengths[i];
        } else {
            lengths.push_back(1);
            prev = *it;
            ++i;
        }
    }

    return lengths;
}

double NumPermsWithRep(const std::vector<int> &v) {
    std::vector<int> myLens = rleCpp(v);
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());

    const int myMax = myLens[0];
    const int numUni = myLens.size();
    double result = 1;

    for (int i = v.size(); i > myMax; --i) {
        result *= i;
    }

    if (numUni > 1) {
        double div = 1.0;

        for (int i = 1; i < numUni; ++i) {
            for (int j = 2; j <= myLens[i]; ++j) {
                div *= j;
            }
        }

        result /= div;
    }

    return result;
}

double NumPermsNoRep(int n, int k) {
    double result = 1.0;

    for (double i = n, m = n - k; i > m; --i) {
        result *= i;
    }

    return result;
}

// The algorithm below is credited to Randy Lai, author of arrangements
// and iterpc. It is much faster than the original naive approach whereby
// we create all combinations of the multiset, thensubsequently count the
// number of permutations of each of those combinations.
double MultisetPermRowNum(int n, int m, const std::vector<int> &Reps) {

    if (n < 2 || m < 1)
        return 1.0;

    const int sumFreqs = std::accumulate(Reps.cbegin(), Reps.cend(), 0);

    if (m == sumFreqs) {
        std::vector<int> freqs(sumFreqs);

        for (int i = 0, k = 0; i < static_cast<int>(Reps.size()); ++i) {
            for (int j = 0; j < Reps[i]; ++j, ++k) {
                freqs[k] = i;
            }
        }

        return NumPermsWithRep(freqs);
    }

    if (m > sumFreqs)
        return 0.0;

    const int n1 = n - 1;
    const int maxFreq = *std::max_element(Reps.cbegin(), Reps.cend());
    const int myMax = (m < maxFreq) ? (m + 1) : (maxFreq + 1);

    // factorial(171)
    // [1] Inf
    // factorial(170)
    // [1] 7.257416e+306
    if (myMax > 170 || m > 170) {
        mpz_t result;
        mpz_init(result);
        MultisetPermRowNumGmp(result, n, m, Reps);

        const double dblRes = (mpz_cmp_d(result, Significand53) > 0) ?
                               std::numeric_limits<double>::infinity() :
                               mpz_get_d(result);

        mpz_clear(result);
        return dblRes;
    }

    std::vector<int> seqR(m);
    std::iota(seqR.begin(), seqR.end(), 1);
    const double prodR = std::accumulate(seqR.cbegin(), seqR.cend(),
                                         1.0, std::multiplies<double>());

    // Create seqeunce from 1 to myMax, then add another
    // 1 at the front... equivalent to c(1, 1:myMax)
    std::vector<double> cumProd(myMax), resV(m + 1, 0.0);
    std::iota(cumProd.begin(), cumProd.end(), 1);

    cumProd.insert(cumProd.begin(), 1);
    std::partial_sum(cumProd.begin(), cumProd.end(),
                     cumProd.begin(), std::multiplies<double>());

    int myMin = std::min(m, Reps[0]);

    for (int i = 0; i <= myMin; ++i) {
        resV[i] = prodR / cumProd[i];
    }

    for (int i = 1; i < n1; ++i) {
        for (int j = m; j > 0; --j) {
            int myMin = std::min(j, Reps[i]);
            double numPerms = 0;

            for (int k = 0; k <= myMin; ++k) {
                numPerms += resV[j - k] / cumProd[k];
            }

            resV[j] = numPerms;
        }
    }

    myMin = std::min(m, Reps[n1]);
    double numPerms = 0.0;

    for (int i = 0; i <= myMin; ++i) {
        numPerms += resV[m - i] / cumProd[i];
    }

    return numPerms;
}

std::vector<int> nonZeroVec(const std::vector<int> &v) {
    std::vector<int> nonZero;

    for (std::size_t i = 0; i < v.size(); i++) {
        if (v[i] > 0) {
            nonZero.push_back(v[i]);
        }
    }

    return nonZero;
}
