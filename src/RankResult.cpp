#include "Permutations/RankPermutation.h"
#include "Combinations/RankCombination.h"
#include "RankResult.h"

rankResultPtr GetRankResultFunc(bool IsComb, bool IsMult,
                               bool IsRep, bool IsGmp) {
    if (IsComb) {
        return GetRankCombFunc(IsMult, IsRep, IsGmp);
    } else {
        return GetRankPermFunc(IsMult, IsRep, IsGmp);
    }
}
