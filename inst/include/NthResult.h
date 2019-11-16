#ifndef NTH_RESULT_H
#define NTH_RESULT_H

#include <Rcpp.h>
#include <gmp.h>

using nthResultPtr = std::vector<int> (*const)(int n, int r, double dblIdx,
                                       mpz_t mpzIdx, const std::vector<int> &Reps);

Rcpp::XPtr<nthResultPtr> putNthResPtrInXPtr(bool IsComb, bool IsMult,
                                            bool IsRep, bool IsGmp);

std::vector<int> nthComb(int n, int r, double dblIdx, 
                         mpz_t mpzIdx, const std::vector<int> &Reps);

std::vector<int> nthCombGmp(int n, int r, double dblIdx, 
                            mpz_t mpzIdx, const std::vector<int> &Reps);

#endif
