#ifndef NTH_RESULT_H
#define NTH_RESULT_H

#include "Permutations/NthPermutation.h"
#include "Combinations/NthCombination.h"

using nthResultPtr = std::vector<int> (*const)(int n, int r,
                                       double dblIdx, mpz_t mpzIdx,
                                       const std::vector<int> &Reps);

nthResultPtr GetNthResultFunc(bool IsComb, bool IsMult,
                              bool IsRep, bool IsGmp) {
    if (IsComb) {
        return GetNthCombFunc(IsMult, IsRep, IsGmp);
    } else {
        return GetNthPermFunc(IsMult, IsRep, IsGmp);
    }
}

#endif