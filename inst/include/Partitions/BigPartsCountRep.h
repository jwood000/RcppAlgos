#ifndef BIG_PARTS_COUNT_REP_H
#define BIG_PARTS_COUNT_REP_H

#include <gmp.h>

void CountPartsRepLenCap(mpz_t res, int n, int m, int cap);
void CountPartsRepLen(mpz_t res, int n, int m);
void CountPartsRep(mpz_t res, int n);

#endif
