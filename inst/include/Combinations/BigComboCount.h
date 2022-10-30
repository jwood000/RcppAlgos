#pragma once

#include <gmpxx.h>
#include <deque>

void nChooseKGmp(mpz_class &result, int n, int m);
void NumCombsWithRepGmp(mpz_class &result, int n, int m);
void MultisetCombRowNumGmp(mpz_class &result, int n, int m,
                           const std::deque<int> &Reps);
