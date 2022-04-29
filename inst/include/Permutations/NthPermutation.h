#ifndef NTH_PERMUTATION_H
#define NTH_PERMUTATION_H

#include <vector>
#include <gmp.h>

using nthPermPtr = std::vector<int> (*const)(int n, int m, double dblIdx,
                                     mpz_t mpzIdx, const std::vector<int> &Reps);

nthPermPtr GetNthPermFunc(bool IsMult, bool IsRep, bool IsGmp);

void TopOffPerm(std::vector<int> &z, const std::vector<int> &myReps,
                int n, int m, bool IsRep, bool IsMult);

#endif
