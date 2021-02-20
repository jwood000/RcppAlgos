#ifndef BIG_PERMUTATION_COUNT_H
#define BIG_PERMUTATION_COUNT_H

#include <gmp.h>
#include <vector>

void NumPermsWithRepGmp(mpz_t result, const std::vector<int> &v);
void NumPermsNoRepGmp(mpz_t result, int n, int k);
void MultisetPermRowNumGmp(mpz_t result, int n, int r,
                           const std::vector<int> &myReps);

#endif
