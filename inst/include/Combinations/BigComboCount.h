#ifndef BIG_COMINATION_COUNT_H
#define BIG_COMINATION_COUNT_H

#include <gmp.h>
#include <deque>

void nChooseKGmp(mpz_t result, int n, int k);
void NumCombsWithRepGmp(mpz_t result, int n, int r);
void MultisetCombRowNumGmp(mpz_t result, int n, int r,
                           const std::deque<int> &Reps);

#endif
