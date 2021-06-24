#ifndef PARTITIONS_COUNT_DISTINCT_H
#define PARTITIONS_COUNT_DISTINCT_H

#include <vector>
#include <gmp.h>

double CountPartDistinctLenCap(int n, int m, int cap);
double CountPartDistinctLen(int n, int m);
double CountPartDistinct(int n);
double CountPartPermDistinct(const std::vector<int> &z,
                             int tar, int width, bool includeZero);
double CountPartPermDistinctCap(const std::vector<int> &z, int cap,
                                int tar, int width, bool includeZero);

void CountPartDistinctLen(mpz_t res, int n, int m);
void CountPartDistinct(mpz_t res, int n);


#endif
