#ifndef BIG_PERMUTATION_COUNT_H
#define BIG_PERMUTATION_COUNT_H

#include <vector>
#include <gmp.h>

void NumPermsWithRepGmp(mpz_t result, const std::vector<int> &v);
void NumPermsNoRepGmp(mpz_t result, int n, int m);
void MultisetPermRowNumGmp(mpz_t result, int n, int m,
                           const std::vector<int> &myReps);

#endif
