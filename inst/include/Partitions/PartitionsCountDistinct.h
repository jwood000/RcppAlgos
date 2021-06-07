#ifndef PARTITIONS_COUNT_DISTINCT_H
#define PARTITIONS_COUNT_DISTINCT_H

#include <vector>

double CountPartDistinctLenCap(int n, int m, int cap);
double CountPartDistinctLen(int n, int m);
double CountPartDistinct(int n);
double CountPartPermDistinct(const std::vector<int> &z,
                             int tar, int width, bool includeZero);

#endif
