#pragma once

#include "CppConvert/GmpxxCopy.h"
#include <vector>

void CountPartsDistinctLenCap(mpz_class &res, std::vector<mpz_class> &p1,
                              std::vector<mpz_class> &p2, int n, int m,
                              int cap, int strtLen);
void CountPartsDistinctLen(mpz_class &res, std::vector<mpz_class> &p1,
                           std::vector<mpz_class> &p2, int n, int m,
                           int cap, int strtLen);
void CountPartsDistinct(mpz_class &res, int n, int m, int cap, int strtLen);
void CountPartsDistinctMultiZero(mpz_class &res, std::vector<mpz_class> &p1,
                                 std::vector<mpz_class> &p2, int n, int m,
                                 int cap, int strtLen);
void CountPartsDistinctCapMZ(mpz_class &res, std::vector<mpz_class> &p1,
                             std::vector<mpz_class> &p2, int n, int m,
                             int cap, int strtLen);
