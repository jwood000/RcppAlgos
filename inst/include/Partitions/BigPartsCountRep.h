#pragma once

#include <gmpxx.h>
#include <vector>

void CountPartsRepLenCap(mpz_class &res, std::vector<mpz_class> &p1,
                         std::vector<mpz_class> &p2, int n, int m,
                         int cap, int strtLen = 0);
void CountPartsRepLen(mpz_class &res, std::vector<mpz_class> &p1,
                      std::vector<mpz_class> &p2, int n, int m,
                      int cap = 0, int strtLen = 0);
void CountPartsRep(mpz_class &res, int n, int m,
                   int cap = 0, int strtLen = 0);
void CountCompsRepLen(mpz_class &res, int n, int m,
                      int cap = 0, int strtLen = 0);
void CountCompsRepZNotWk(mpz_class &res, int n, int m,
                         int cap = 0, int strtLen = 0);
