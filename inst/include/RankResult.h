#pragma once

#include <vector>
#include "CppConvert/GmpxxCopy.h"

using rankResultPtr = void (*const)(std::vector<int>::iterator iter, int n,
                            int m, double &dblIdx, mpz_class &mpzIdx,
                            const std::vector<int> &Reps);

rankResultPtr GetRankResultFunc(bool IsComb, bool IsMult,
                                bool IsRep, bool IsGmp);
