#ifndef SET_UP_UTILS_H
#define SET_UP_UTILS_H

#include "Combinations/NthCombination.h"
#include "Permutations/NthPermutation.h"
#include "CleanConvert.h"

void SetType(VecType &myType, SEXP Rv);
void SetFactorClass(SEXP res, SEXP Rv);

void SetValues(VecType &myType, std::vector<int> &Reps,
               std::vector<int> &freqs, std::vector<int> &vInt,
               std::vector<double> &vNum, SEXP Rv, SEXP RFreqs,
               SEXP Rm, int &n, int &m, bool &IsMult,
               bool &IsRep, bool IsConstrained = false);

void SetThreads(bool &Parallel, int maxThreads, int nRows,
                VecType myType, int &nThreads, SEXP RNumThreads, int limit);

void SetNumResults(bool IsGmp, bool bLower, bool bUpper, bool IsGenCnstrd,
                   mpz_t *const upperMpz, mpz_t *const lowerMpz, double lower,
                   double upper, double computedRows, mpz_t &computedRowsMpz,
                   int &nRows, double &userNumRows);

void SetBounds(SEXP Rlow, SEXP Rhigh, bool IsGmp, bool &bLower,
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

void SetRandomSample(SEXP RindexVec, SEXP RNumSamp, int &sampSize,
                     bool IsGmp, double computedRows,
                     std::vector<double> &mySample,
                     SEXP baseSample, SEXP rho, bool bSubOne = true);

void SetRandomSampleMpz(const SEXP &RindexVec, const SEXP &RmySeed,
                        int sampSize, bool IsGmp, mpz_t &computedRowsMpz,
                        mpz_t *const myVec, bool bSubOne = true);

void SetSampleNames(SEXP objRcpp, bool IsGmp, int sampSize,
                    const std::vector<double> &mySample,
                    mpz_t *const myBigSamp, bool bAddOne = true);

#endif
