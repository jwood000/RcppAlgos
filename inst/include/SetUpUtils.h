#ifndef SET_UP_UTILS_H
#define SET_UP_UTILS_H

#include "Combinations/NthCombination.h"
#include "Permutations/NthPermutation.h"
#include "CleanConvert.h"
#include <gmp.h>

void SetType(VecType &myType, const SEXP &Rv);
void SetFactorClass(SEXP res, const SEXP &Rv);
void SetValues(VecType &myType, std::vector<int> &vInt,
               std::vector<double> &vNum, int &n, const SEXP &Rv);

void SetFreqsAndM(SEXP RFreqs, bool &IsMult, std::vector<int> &Reps, bool &IsRep,
                  std::vector<int> &freqs, const SEXP &Rm, int n, int &m);

void SetThreads(bool &Parallel, int maxThreads, int nRows,
                VecType myType, int &nThreads, SEXP RNumThreads, int limit);

void SetNumResults(bool IsGmp, bool bLower, bool bUpper, bool IsGenCnstrd,
                   mpz_t *const upperMpz, mpz_t *const lowerMpz, double lower,
                   double upper, double computedRows, mpz_t &computedRowsMpz,
                   int &nRows, double &userNumRows);

void SetBounds(const SEXP &Rlow, const SEXP &Rhigh, bool IsGmp, bool &bLower,
               bool &bUpper, double &lower, double &upper, mpz_t *const lowerMpz,
               mpz_t *const upperMpz, mpz_t computedRowMpz, double computedRows);

void SetStartZ(const std::vector<int> &myReps,
               const std::vector<int> &freqs,
               std::vector<int> &z, bool IsComb, int n,
               int m, double lower, mpz_t lowerMpz,
               bool IsRep, bool IsMult, bool IsGmp);

void PermuteSpecific(int &phaseOne, bool &generalRet, int n, int m,
                     int nRows, bool IsMult, bool IsCharacter,
                     bool IsComb, bool bLower, bool IsRep);
    
#endif
