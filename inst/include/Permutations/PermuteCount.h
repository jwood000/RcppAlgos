#pragma once

#include <vector>

std::vector<int> rleCpp(const std::vector<int> &x, int first_idx = 0);
double NumPermsWithRep(const std::vector<int> &v, bool includeZero = true);
double NumPermsNoRep(int n, int m);
double MultisetPermRowNum(int n, int m, const std::vector<int> &Reps);
std::vector<int> nonZeroVec(const std::vector<int> &v);
