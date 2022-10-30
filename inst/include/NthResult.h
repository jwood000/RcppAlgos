#pragma once

#include <vector>
#include <gmpxx.h>

using nthResultPtr = std::vector<int> (*const)(int n, int m,
                                       double dblIdx, const mpz_class &mpzIdx,
                                       const std::vector<int> &Reps);

nthResultPtr GetNthResultFunc(bool IsComb, bool IsMult,
                              bool IsRep, bool IsGmp);

void SetNextIter(const std::vector<int> &myReps, std::vector<int> &z,
                 const nthResultPtr nthRes, double &lower, mpz_class &lowerMpz,
                 int stepSize, int n, int m, bool IsGmp, bool IsComb,
                 bool IsRep, bool IsMult);
