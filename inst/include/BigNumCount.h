#ifndef BIG_NUM_COUNT_H
#define BIG_NUM_COUNT_H

#include <gmp.h>
#include <vector>

void NumPermsWithRepGmp(mpz_t result, const std::vector<int> &v);
void NumPermsNoRepGmp(mpz_t result, int n, int k);
void nChooseKGmp(mpz_t result, int n, int k);
void NumCombsWithRepGmp(mpz_t result, int n, int r);
void MultisetCombRowNumGmp(mpz_t result, int n, int r, const std::vector<int> &Reps);
void MultisetPermRowNumGmp(mpz_t result, int n, int r, const std::vector<int> &myReps);

void GetComputedRowMpz(mpz_t computedRowMpz, bool IsMultiset, bool IsComb, bool IsRep,
                       int n, int m, const SEXP &Rm, const std::vector<int> &freqs, 
                       const std::vector<int> &myReps);

#endif
