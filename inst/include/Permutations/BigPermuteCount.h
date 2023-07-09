#pragma once

#include <vector>
#include <gmpxx.h>

void NumPermsWithRepGmp(mpz_class &result, const std::vector<int> &v);
void NumPermsNoRepGmp(mpz_class &result, int n, int m);
void MultisetPermRowNumGmp(mpz_class &result, int n, int m,
                           const std::vector<int> &myReps);
