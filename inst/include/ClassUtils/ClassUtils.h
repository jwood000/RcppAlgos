#pragma once

#include "Constraints/ConstraintsTypes.h"
#include "CppConvert.h"

void SetIndexVec(SEXP RindexVec, std::vector<double> &mySample,
                 std::size_t &sampSize, bool IsGmp, double computedRows);

void SetIndexVecMpz(SEXP RindexVec, std::vector<mpz_class> &myVec,
                    std::size_t sampSize, mpz_class computedRowsMpz);

void increment(bool IsGmp, mpz_class &mpzIndex, double &dblIndex);
void increment(bool IsGmp, mpz_class &mpzIndex, double &dblIndex, int nRows);
void decrement(bool IsGmp, mpz_class &mpzIndex, double &dblIndex);
void decrement(bool IsGmp, mpz_class &mpzIndex, double &dblIndex, int nRows);

bool CheckEqSi(bool IsGmp, const mpz_class &mpzIndex, double dblIndex, int si);
bool CheckIndLT(bool IsGmp, const mpz_class &mpzIndex, double dblIndex,
                const mpz_class &computedRowsMpz, double computedRows,
                bool eq = false);
bool CheckEqInd(bool IsGmp, const mpz_class &mpzIndex, double dblIndex,
                const mpz_class &computedRowsMpz, double computedRows);
bool CheckIndGrT(bool IsGmp, const mpz_class &mpzIndex, double dblIndex,
                 const mpz_class &computedRowsMpz, double computedRows);
bool CheckGrTSi(bool IsGmp, const mpz_class &mpzIndex, double dblIndex, int si);

void zUpdateIndex(const std::vector<double> &vNum,
                  const std::vector<int> &vInt, std::vector<int> &z,
                  SEXP v, SEXP mat, std::size_t m, std::size_t nRows,
                  bool bAddOne = false);
