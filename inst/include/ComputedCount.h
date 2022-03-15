#ifndef COMPUTED_COUNT_H
#define COMPUTED_COUNT_H

#include "CleanConvert.h"

double GetComputedRows(bool IsMult, bool IsComb, bool IsRep,
                       int n, int m, const SEXP &Rm,
                       const std::vector<int> &freqs,
                       const std::vector<int> &Reps);

void GetComputedRowMpz(mpz_t computedRowMpz, bool IsMult, bool IsComb,
                       bool IsRep, int n, int m, const SEXP &Rm,
                       const std::vector<int> &freqs,
                       const std::vector<int> &myReps);

#endif