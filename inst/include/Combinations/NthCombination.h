#pragma once

#include <vector>
#include "CppConvert/GmpxxCopy.h"

using nthCombPtr = std::vector<int> (*const)(int n, int m, double dblIdx,
                                     const mpz_class &mpzIdx,
                                     const std::vector<int> &Reps);

nthCombPtr GetNthCombFunc(bool IsMult, bool IsRep, bool IsGmp);

std::vector<int> nthComb(int n, int m, double dblIdx,
                         const mpz_class &mpzIdx,
                         const std::vector<int> &Reps);

std::vector<int> nthCombGmp(int n, int m, double dblIdx,
                            const mpz_class &mpzIdx,
                            const std::vector<int> &Reps);
