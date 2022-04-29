#ifndef RANK_PERMUTATION_H
#define RANK_PERMUTATION_H

#include <vector>
#include <gmp.h>

using It = std::vector<int>::iterator;

using rankPermPtr = void (*const)(It iter, int n, int r, double &dblIdx,
                          mpz_t mpzIdx, const std::vector<int> &Reps);

rankPermPtr GetRankPermFunc(bool IsMult, bool IsRep, bool IsGmp);

#endif
