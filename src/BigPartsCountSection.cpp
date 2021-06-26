#include <gmp.h>

void SumSection(mpz_t n, mpz_t res) {
    mpz_t nIter;
    mpz_t sumOne;
    mpz_t sumTwo;

    mpz_init(nIter);
    mpz_init(sumOne);
    mpz_init(sumTwo);

    mpz_div_ui(nIter, n, 3u);

    mpz_sub_ui(sumOne, n, 1u);
    mpz_mul(sumOne, sumOne, nIter);

    mpz_sub_ui(sumTwo, nIter, 1u);
    mpz_mul(sumTwo, sumTwo, nIter);
    mpz_div_2exp(sumTwo, sumTwo, 1);
    mpz_mul_ui(sumTwo, sumTwo, 3u);

    mpz_div_2exp(res, nIter, 1);
    mpz_add(res, res, sumTwo);
    mpz_sub(res, sumOne, res);
    mpz_div_2exp(res, res, 1);

    mpz_clear(nIter);
    mpz_clear(sumOne);
    mpz_clear(sumTwo);
}
