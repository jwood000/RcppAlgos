#ifndef RcppAlgos_NthResult_h
#define RcppAlgos_NthResult_h

#include <CombPermUtility.h>

std::vector<int> nthCombination(int n, int r, double myIndex, bool isRep,
                                bool isMult, std::vector<int> Reps);

std::vector<int> nthPermutation(int n, int r, double myIndex, bool isRep,
                                bool isMult, std::vector<int> Reps);

#endif
