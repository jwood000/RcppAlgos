#ifndef RANK_RESULT_H
#define RANK_RESULT_H

#include <vector>
#include <gmp.h>

using rankResultPtr = void (*const)(std::vector<int>::iterator iter, int n,
                            int m, double &dblIdx, mpz_t mpzIdx,
                            const std::vector<int> &Reps);

rankResultPtr GetRankResultFunc(bool IsComb, bool IsMult,
                                bool IsRep, bool IsGmp);

#endif