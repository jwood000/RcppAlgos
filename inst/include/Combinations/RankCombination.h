#ifndef RANK_COMBINATION_H
#define RANK_COMBINATION_H

#include <vector>
#include <gmp.h>

using It = std::vector<int>::iterator;

using rankCombPtr = void (*const)(It iter, int n, int r, double &dblIdx,
                          mpz_t mpzIdx, const std::vector<int> &Reps);

rankCombPtr GetRankCombFunc(bool IsMult, bool IsRep, bool IsGmp);

#endif
