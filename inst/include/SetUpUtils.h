#pragma once

#include "CppConvert.h"

void SetType(VecType &myType, SEXP Rv);
void SetFactorClass(SEXP res, SEXP Rv);
int GetLength(SEXP Rv, VecType myType);

void SetFreqsAndM(std::vector<int> &Reps,
                  std::vector<int> &freqs, SEXP RFreqs, SEXP Rm,
                  int &n, int &m, bool &IsMult, bool &IsRep);

void SetBasic(SEXP Rv, std::vector<double> &vNum,
              std::vector<int> &vInt, int &n, VecType &myType);

void SetValues(VecType &myType, std::vector<int> &Reps,
               std::vector<int> &freqs, std::vector<int> &vInt,
               std::vector<double> &vNum, SEXP Rv, SEXP RFreqs,
               SEXP Rm, int &n, int &m, bool &IsMult,
               bool &IsRep, bool IsConstrained = false);

void SetThreads(bool &Parallel, int maxThreads, int nRows,
                VecType myType, int &nThreads, SEXP RNumThreads, int limit);

void SetNumResults(bool IsGmp, bool bLower, bool bUpper, bool bSetNum,
                   const mpz_class &upperMpz, const mpz_class &lowerMpz,
                   double lower, double upper, double computedRows,
                   const mpz_class &computedRowsMpz, int &nRows,
                   double &userNumRows);

void SetBounds(SEXP Rlow, SEXP Rhigh, bool IsGmp, bool &bLower,
               bool &bUpper, double &lower, double &upper,
               mpz_class &lowerMpz, mpz_class &upperMpz,
               const mpz_class &computedRowMpz, double computedRows);

void SetStartZ(const std::vector<int> &myReps,
               const std::vector<int> &freqs, std::vector<int> &z,
               bool IsComb, int n, int m, double lower,
               const mpz_class &lowerMpz, bool IsRep,
               bool IsMult, bool IsGmp);

void PermuteSpecific(int &phaseOne, bool &generalRet, int n, int m,
                     int nRows, bool IsMult, bool IsCharacter,
                     bool IsComb, bool bLower, bool IsRep);

void SetRandomSample(SEXP RindexVec, SEXP RNumSamp, int &sampSize,
                     bool IsGmp, double computedRows,
                     std::vector<double> &mySample,
                     SEXP baseSample, SEXP rho);

void SetRandomSampleMpz(SEXP RindexVec, SEXP RmySeed, int sampSize,
                        bool IsGmp, const mpz_class &computedRowsMpz,
                        std::vector<mpz_class> &myVec);

void SetSampleNames(SEXP objRcpp, bool IsGmp, int sampSize,
                    const std::vector<double> &mySample,
                    const std::vector<mpz_class> &myBigSamp, bool IsNamed,
                    SEXP colNames = R_NilValue, int xtraDims = 0);

SEXP GetInt64Vec(const std::vector<std::int64_t> &v);
