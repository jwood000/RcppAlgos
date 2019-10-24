#ifndef NTH_RESULT_H
#define NTH_RESULT_H

#include <Rcpp.h>
#include <gmp.h>

using nthResultPtr = std::vector<int> (*)(int n, int r, double dblIdx,
                                       mpz_t mpzIdx, std::vector<int> Reps);

Rcpp::XPtr<nthResultPtr> putNthResPtrInXPtr(bool IsMultiset, bool IsRep,
                                            bool IsGmp, bool IsComb);

std::vector<int> nthComb(int n, int r, double dblIdx, 
                         mpz_t mpzIdx, std::vector<int> Reps);

std::vector<int> nthCombGmp(int n, int r, double dblIdx, 
                            mpz_t mpzIdx, std::vector<int> Reps);

#endif
