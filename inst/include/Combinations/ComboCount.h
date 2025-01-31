#pragma once

#include <vector>
#include <deque>

double nChooseK(int n, int m);
double NumCombsWithRep(int n, int m);
double MultisetCombRowNumFast(int n,int m, const std::vector<int> &Reps);

// This one isn't as efficient as MultisetCombRowNumFast, however it will
// not produce negative results and is thus used in determining whether
// gmp analogs are necessary
double MultisetCombRowNum(int n, int m, const std::vector<int> &Reps);

void ManageCountsVector(std::vector<int> &Counts, int &n1);
void ManageCountsDeque(std::deque<int> &Counts, int &n1);
