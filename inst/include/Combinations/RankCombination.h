#pragma once

#include <vector>
#include "CppConvert/GmpxxCopy.h"

using rankCombPtr = void (*const)(std::vector<int>::iterator iter, int n,
                          int m, double &dblIdx, mpz_class &mpzIdx,
                          const std::vector<int> &Reps);

rankCombPtr GetRankCombFunc(bool IsMult, bool IsRep, bool IsGmp);
