#ifndef RANK_PERMUTATION_H
#define RANK_PERMUTATION_H

#include <vector>
#include <gmp.h>

using rankPermPtr = void (*const)(std::vector<int>::iterator iter, int n,
                          int m, double &dblIdx, mpz_t mpzIdx,
                          const std::vector<int> &Reps);

rankPermPtr GetRankPermFunc(bool IsMult, bool IsRep, bool IsGmp);

#endif
