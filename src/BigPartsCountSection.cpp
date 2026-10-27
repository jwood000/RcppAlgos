#include "Partitions/BigPartsCountSection.h"
#include <algorithm>

void SumSection(const mpz_class &n, mpz_class &res) {
    mpz_class nIter(n / 3);
    mpz_class sumOne;
    mpz_class sumTwo;

    sumOne = (n - 1) * nIter;
    sumTwo = (nIter - 1) * nIter;
    mpz_div_2exp(sumTwo.get_mpz_t(), sumTwo.get_mpz_t(), 1);
    sumTwo *= 3;

    mpz_div_2exp(res.get_mpz_t(), nIter.get_mpz_t(), 1);
    res += sumTwo;
    res = sumOne - res;
    mpz_div_2exp(res.get_mpz_t(), res.get_mpz_t(), 1);
}

void ResetP1D(std::vector<mpz_class> &p) {
    std::fill(p.begin(), p.end(), 0);
}

void ResetP2D(std::vector<std::vector<mpz_class>> &p2d) {
    for (auto& v : p2d) {
        ResetP1D(v);
    }
}
