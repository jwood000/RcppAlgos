#include "Permutations/RankPermutation.h"
#include "Combinations/RankCombination.h"

using rankResultPtr = void (*const)(std::vector<int>::iterator iter, int n,
                            int m, double &dblIdx, mpz_t mpzIdx,
                            const std::vector<int> &Reps);

rankResultPtr GetRankResultFunc(bool IsComb, bool IsMult,
                               bool IsRep, bool IsGmp) {
    if (IsComb) {
        return GetRankCombFunc(IsMult, IsRep, IsGmp);
    } else {
        return GetRankPermFunc(IsMult, IsRep, IsGmp);
    }
}
