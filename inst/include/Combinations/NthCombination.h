#ifndef NTH_COMBINATION_H
#define NTH_COMBINATION_H

#include <vector>
#include <gmp.h>

using nthCombPtr = std::vector<int> (*const)(int n, int r, double dblIdx,
                                     mpz_t mpzIdx, const std::vector<int> &Reps);

nthCombPtr GetNthCombFunc(bool IsMult, bool IsRep, bool IsGmp);

std::vector<int> nthComb(int n, int r, double dblIdx,
                         mpz_t mpzIdx, const std::vector<int> &Reps);

std::vector<int> nthCombGmp(int n, int r, double dblIdx,
                            mpz_t mpzIdx, const std::vector<int> &Reps);

#endif
