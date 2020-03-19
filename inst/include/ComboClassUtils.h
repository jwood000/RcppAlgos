#ifndef COMBO_CLASS_UTILS_H
#define COMBO_CLASS_UTILS_H

#include "CleanConvert.h"

void TopOffPartialPerm(std::vector<int> &z, const std::vector<int> &myReps,
                       int n, int m, bool IsComb, bool IsRep, bool IsMult);

void SetIndexVec(SEXP RindexVec, std::vector<double> &mySample,
                 std::size_t &sampSize, bool IsGmp, double computedRows);

void SetIndexVecMpz(SEXP RindexVec, mpz_t *myVec,
                    std::size_t sampSize, mpz_t &computedRowsMpz);

void increment(bool IsGmp, mpz_t &mpzIndex, double &dblIndex);
void increment(bool IsGmp, mpz_t &mpzIndex, double &dblIndex, int nRows);
void decrement(bool IsGmp, mpz_t &mpzIndex, double &dblIndex);
void decrement(bool IsGmp, mpz_t &mpzIndex, double &dblIndex, int nRows);

bool CheckEqSi(bool IsGmp, mpz_t &mpzIndex, double &dblIndex, int si);
bool CheckIndLT(bool IsGmp, mpz_t &mpzIndex, double &dblIndex,
                mpz_t computedRowsMpz, double computedRows, bool eq = false);
bool CheckEqInd(bool IsGmp, mpz_t &mpzIndex, double &dblIndex,
                mpz_t computedRowsMpz, double computedRows);
bool CheckIndGrT(bool IsGmp, mpz_t &mpzIndex, double &dblIndex,
                 mpz_t computedRowsMpz, double computedRows);
bool CheckGrTSi(bool IsGmp, mpz_t &mpzIndex, double &dblIndex, int si);

#endif
