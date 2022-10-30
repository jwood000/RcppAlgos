#pragma once

#include <vector>
#include <gmpxx.h>

using rankPermPtr = void (*const)(std::vector<int>::iterator iter, int n,
                          int m, double &dblIdx, mpz_class &mpzIdx,
                          const std::vector<int> &Reps);

rankPermPtr GetRankPermFunc(bool IsMult, bool IsRep, bool IsGmp);
