#ifndef BIG_PARTS_COUNT_DISTINCT_H
#define BIG_PARTS_COUNT_DISTINCT_H

#include <gmp.h>

void CountPartsDistinctLenCap(mpz_t res, mpz_t* p1, mpz_t* p2,
                              int n, int m, int cap, int strtLen);
void CountPartsDistinctLen(mpz_t res, mpz_t* p1, mpz_t* p2,
                           int n, int m, int cap, int strtLen);
void CountPartsDistinct(mpz_t res, int n, int m, int cap, int strtLen);
void CountPartsDistinctMultiZero(mpz_t res, mpz_t* p1, mpz_t* p2,
                                 int n, int m, int cap, int strtLen);
void CountPartsDistinctCapMZ(mpz_t res, mpz_t* p1, mpz_t* p2,
                             int n, int m, int cap, int strtLen);

#endif
