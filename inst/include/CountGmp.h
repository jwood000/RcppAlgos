#ifndef COUNT_GMP_H
#define COUNT_GMP_H

#include <Rcpp.h>
#include <gmp.h>

void NumPermsWithRepGmp(mpz_t result, std::vector<int> &v);
void NumPermsNoRepGmp(mpz_t result, int n, int k);
void nChooseKGmp(mpz_t result, int n, int k);
void NumCombsWithRepGmp(mpz_t result, int n, int r);
void MultisetCombRowNumGmp(mpz_t result, int n, int r, std::vector<int> &Reps);
void MultisetPermRowNumGmp(mpz_t result, int n, int r, std::vector<int> &myReps);

#endif
