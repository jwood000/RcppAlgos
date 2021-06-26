#ifndef BIG_PARTS_COUNT_DISTINCT_H
#define BIG_PARTS_COUNT_DISTINCT_H

#include <gmp.h>

void CountPartsDistinctLenCap(mpz_t res, int n, int m, int cap);
void CountPartsDistinctLen(mpz_t res, int n, int m);
void CountPartsDistinct(mpz_t res, int n);
void CountPartsDistinctMultiZero(mpz_t res, int target, int m, int strtLen);

#endif
