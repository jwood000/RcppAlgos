#ifndef PARTITIONS_COUNT_REP_H
#define PARTITIONS_COUNT_REP_H

#include <gmp.h>

double CountPartRepLenCap(int n, int m, int cap);
double CountPartRepLen(int n, int m);
double CountPartRep(int n);
double CountPartPermRep(int target, int m, bool includeZero);

void CountPartRepLen(mpz_t res, int n, int m);
void CountPartRep(mpz_t res, int n);

#endif
