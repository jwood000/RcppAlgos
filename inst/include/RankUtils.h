#ifndef RANK_UTILS_H
#define RANK_UTILS_H

#include "SetUpUtils.h"

void SetUpRank(SEXP RIdx, SEXP Rv, SEXP RisRep, SEXP RFreqs, SEXP Rm,
               std::vector<int> &idx, std::vector<int> &freqs,
               std::vector<int> &myReps, VecType &myType, int &n, int &m,
               bool IsComb, bool &IsMult, bool &IsRep);

SEXP MpzReturn(mpz_t *myVec, int numResults);

#endif