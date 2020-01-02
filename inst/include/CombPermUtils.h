#ifndef COMB_PERM_UTILS_H
#define COMB_PERM_UTILS_H

#include "CleanConvert.h"

void SetType(VecType &myType, const SEXP &Rv);
std::vector<int> zUpdateIndex(SEXP x, SEXP y, int m);

void SetFactorClass(Rcpp::IntegerMatrix &matInt, const SEXP &Rv);

void SetValues(VecType &myType, std::vector<int> &vInt,
               std::vector<double> &vNum, int &n, const SEXP &Rv);

void SetFreqsAndM(SEXP RFreqs, bool &IsMultiset, std::vector<int> &Reps, bool &IsRepetition,
                  int &lenFreqs, std::vector<int> &freqsExpanded, const SEXP &Rm, int n, int &m);

void SetThreads(bool &Parallel, int maxThreads, int nRows,
                VecType myType, int &nThreads, SEXP RNumThreads, int limit);

SEXP CopyRv(const SEXP &Rv, const std::vector<int> &vInt,
            const std::vector<double> &vNum, VecType myType, bool IsFactor = false);

void SetRandomSample(SEXP RindexVec, SEXP RNumSamp, std::size_t &sampSize,
                     bool IsGmp, double computedRows, std::vector<double> &mySample,
                     Rcpp::Function baseSample);

std::vector<int> rleCpp(const std::vector<int> &x);
double NumPermsWithRep(const std::vector<int> &v);
double NumPermsNoRep(int n, int k);
double nChooseK(int n, int k);
double NumCombsWithRep(int n, int r);
double MultisetCombRowNumFast(int n,int r, const std::vector<int> &Reps);
double MultisetPermRowNum(int n, int r, const std::vector<int> &Reps);

// This one isn't as efficient as MultisetCombRowNumFast, however it will
// not produce negative results and is thus used in determining whether
// gmp analogs are necessary
double MultisetCombRowNum(int n, int r, const std::vector<int> &Reps);

double GetComputedRows(bool IsMultiset, bool IsComb, bool IsRep, int n, int &m, SEXP Rm,
                       int lenFreqs, std::vector<int> &freqs, std::vector<int> &Reps);
    
#endif
