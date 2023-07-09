#pragma once

#include <gmpxx.h>
#include <vector>

void CountPartsRepLenCap(mpz_class &res, std::vector<mpz_class> &p1,
                         std::vector<mpz_class> &p2, int n, int m,
                         int cap, int strtLen);
void CountPartsRepLen(mpz_class &res, std::vector<mpz_class> &p1,
                      std::vector<mpz_class> &p2, int n, int m,
                      int cap, int strtLen);
void CountPartsRep(mpz_class &res, int n, int m, int cap, int strtLen);
void CountCompsRepLen(mpz_class &res, int n, int m, int cap, int strtLen);
void CountCompsRepZero(mpz_class &res, int n, int m, int cap, int strtLen);
