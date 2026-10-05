#pragma once

#include <vector>

// v must be sorted
double NumPermsWithRep(const std::vector<int> &v);
double NumPermsNoRep(int n, int m);
double MultisetPermRowNum(int n, int m, const std::vector<int> &Reps);
std::vector<int> nonZeroVec(const std::vector<int> &v);
