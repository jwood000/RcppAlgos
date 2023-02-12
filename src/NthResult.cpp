#include "CppConvert/GmpxxCopy.h"
#include <vector>

#include "Permutations/NthPermutation.h"
#include "Combinations/NthCombination.h"

using nthResultPtr = std::vector<int> (*const)(int n, int m,
                                       double dblIdx, const mpz_class &mpzIdx,
                                       const std::vector<int> &Reps);

nthResultPtr GetNthResultFunc(bool IsComb, bool IsMult,
                              bool IsRep, bool IsGmp) {
    if (IsComb) {
        return GetNthCombFunc(IsMult, IsRep, IsGmp);
    } else {
        return GetNthPermFunc(IsMult, IsRep, IsGmp);
    }
}

void SetNextIter(const std::vector<int> &myReps, std::vector<int> &z,
                 const nthResultPtr nthResFun, double &lower,
                 mpz_class &lowerMpz, int stepSize, int n, int m, bool IsGmp,
                 bool IsComb, bool IsRep, bool IsMult) {

    if (IsGmp) {
        lowerMpz += stepSize;
    } else {
        lower += stepSize;
    }

    z = nthResFun(n, m, lower, lowerMpz, myReps);

    if (!IsComb) {
        TopOffPerm(z, myReps, n, m, IsRep, IsMult);
    }
}
