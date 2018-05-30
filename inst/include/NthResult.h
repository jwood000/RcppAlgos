#ifndef NTH_RESULT_H
#define NTH_RESULT_H

#include <vector>

std::vector<int> nthCombination(int n, int r, double myIndex, bool isRep,
                                bool isMult, std::vector<int> Reps);

std::vector<int> nthPermutation(int n, int r, double myIndex, bool isRep,
                                bool isMult, std::vector<int> Reps);

#endif
