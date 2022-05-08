#ifndef BIG_COMINATION_COUNT_H
#define BIG_COMINATION_COUNT_H

#include <gmp.h>
#include <deque>

void nChooseKGmp(mpz_t result, int n, int m);
void NumCombsWithRepGmp(mpz_t result, int n, int m);
void MultisetCombRowNumGmp(mpz_t result, int n, int m,
                           const std::deque<int> &Reps);

#endif
