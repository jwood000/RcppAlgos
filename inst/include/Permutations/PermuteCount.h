#pragma once

#include <vector>

// `v` must be sorted in non-decreasing order. When includeZero is false,
// only leading zeros are excluded from the count.
double NumPermsWithRep(const std::vector<int> &v, bool includeZero = true);
double NumPermsNoRep(int n, int m);
double MultisetPermRowNum(int n, int m, const std::vector<int> &Reps);
std::vector<int> nonZeroVec(const std::vector<int> &v);
