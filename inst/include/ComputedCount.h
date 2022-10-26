#pragma once

#include "CppConvert.h"

double GetComputedRows(bool IsMult, bool IsComb, bool IsRep,
                       int n, int m, const SEXP &Rm,
                       const std::vector<int> &freqs,
                       const std::vector<int> &Reps);

void GetComputedRowMpz(mpz_class &computedRowMpz, bool IsMult,
                       bool IsComb, bool IsRep, int n, int m,
                       const SEXP &Rm, const std::vector<int> &freqs,
                       const std::vector<int> &myReps);
