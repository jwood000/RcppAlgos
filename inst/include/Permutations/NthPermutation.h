#pragma once

#include <vector>
#include "CppConvert/GmpxxCopy.h"

using nthPermPtr = std::vector<int> (*const)(int n, int m, double dblIdx,
                                     const mpz_class &mpzIdx,
                                     const std::vector<int> &Reps);

nthPermPtr GetNthPermFunc(bool IsMult, bool IsRep, bool IsGmp);

void TopOffPerm(std::vector<int> &z, const std::vector<int> &myReps,
                int n, int m, bool IsRep, bool IsMult);
