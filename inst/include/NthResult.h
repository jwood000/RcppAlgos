#ifndef NTH_RESULT_H
#define NTH_RESULT_H

#include <vector>
#include <gmp.h>

using nthResultPtr = std::vector<int> (*const)(int n, int m,
                                       double dblIdx, mpz_t mpzIdx,
                                       const std::vector<int> &Reps);

nthResultPtr GetNthResultFunc(bool IsComb, bool IsMult,
                              bool IsRep, bool IsGmp);

void SetNextIter(const std::vector<int> &myReps, std::vector<int> &z,
                 const nthResultPtr nthRes, double &lower, mpz_t lowerMpz,
                 int stepSize, int n, int m, bool IsGmp, bool IsComb,
                 bool IsRep, bool IsMult);

#endif