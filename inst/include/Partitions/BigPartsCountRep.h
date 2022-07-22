#ifndef BIG_PARTS_COUNT_REP_H
#define BIG_PARTS_COUNT_REP_H

#include <gmp.h>

void CountPartsRepLenCap(mpz_t res, mpz_t* p1, mpz_t* p2,
                         int n, int m, int cap, int strtLen);
void CountPartsRepLen(mpz_t res, mpz_t* p1, mpz_t* p2,
                      int n, int m, int cap, int strtLen);
void CountPartsRep(mpz_t res, int n, int m, int cap, int strtLen);
void CountCompsRepLen(mpz_t res, int n, int m, int cap, int strtLen);
void CountCompsRepZero(mpz_t res, int n, int m, int cap, int strtLen);

#endif
