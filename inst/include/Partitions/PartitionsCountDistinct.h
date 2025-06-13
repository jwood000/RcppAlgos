#pragma once

#include <vector>

double CountPartsDistinctLenCap(int n, int m, int cap, int strtLen = 0);
double CountPartsDistinctLen(int n, int m, int cap = 0, int strtLen = 0);
double CountPartsDistinct(int n, int m, int cap = 0, int strtLen = 0);
double CountPartsDistinctMultiZero(int n, int m, int cap, int strtLen);
double CountPartsDistinctCapMZ(int n, int m, int cap, int strtLen);
double CountPartsPermDistinctCap(int n, int m, int cap, int strtLen = 0);
double CountPartsPermDistinctCapMZ(int n, int m, int cap, int strtLen);
double CountCompsDistinctLen(int n, int m, int cap = 0, int strtLen = 0);
double CountCompsDistinctMultiZero(int n, int m, int cap, int strtLen);
double CountCompsDistinctMZWeak(int n, int m, int cap, int strtLen);
