#ifndef RANK_COMBINATION_H
#define RANK_COMBINATION_H

#include <vector>
#include <gmp.h>

using rankCombPtr = void (*const)(std::vector<int>::iterator iter, int n,
                          int m, double &dblIdx, mpz_t mpzIdx,
                          const std::vector<int> &Reps);

rankCombPtr GetRankCombFunc(bool IsMult, bool IsRep, bool IsGmp);

#endif
