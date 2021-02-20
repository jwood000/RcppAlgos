#ifndef PERMUTE_COUNT_H
#define PERMUTE_COUNT_H

#include <vector>

std::vector<int> rleCpp(const std::vector<int> &x);
double NumPermsWithRep(const std::vector<int> &v);
double NumPermsNoRep(int n, int k);
double MultisetPermRowNum(int n, int r, const std::vector<int> &Reps);

#endif
