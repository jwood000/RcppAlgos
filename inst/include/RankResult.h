#ifndef RANK_RESULT_H
#define RANK_RESULT_H

#include <vector>
#include <gmp.h>

using It = std::vector<int>::iterator;

using rankResultPtr = void (*const)(It iter, int n, int r, double &dblIdx,
                            mpz_t mpzIdx, const std::vector<int> &Reps);

rankResultPtr GetRankResultFunc(bool IsComb, bool IsMult,
                                bool IsRep, bool IsGmp);

#endif