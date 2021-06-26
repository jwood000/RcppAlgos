#ifndef PARTITIONS_COUNT_DISTINCT_H
#define PARTITIONS_COUNT_DISTINCT_H

#include <vector>

double CountPartsDistinctLenCap(int n, int m, int cap);
double CountPartsDistinctLen(int n, int m);
double CountPartsDistinct(int n);
double CountPartsDistinctMultiZero(int target, int m, int strtLen);

double CountPartsPermDistinct(const std::vector<int> &z,
                             int tar, int width, bool includeZero);
double CountPartsPermDistinctCap(const std::vector<int> &z, int cap,
                                int tar, int width, bool includeZero);

#endif
