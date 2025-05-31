#pragma once

#include <vector>

double CountPartsDistinctLenCap(int n, int m, int cap, int strtLen = 0);
double CountPartsDistinctLen(int n, int m, int cap = 0, int strtLen = 0);
double CountPartsDistinct(int n, int m, int cap = 0, int strtLen = 0);
double CountPartsDistinctMultiZero(int n, int m, int cap, int strtLen);
double CountPartsDistinctCapMZ(int n, int m, int cap, int strtLen);

double CountPartsPermDistinct(const std::vector<int> &z,
                             int tar, int width, bool includeZero = false);
double CountPartsPermDistinctCap(const std::vector<int> &z, int cap,
                                int tar, int width, bool includeZero);
double CountCompsDistinctLen(int n, int m, int cap = 0, int strtLen = 0);
double CountCompsDistinctMultiZero(int n, int m, int cap, int strtLen);
double CountCompsDistinctMZWeak(int n, int m, int cap, int strtLen);
