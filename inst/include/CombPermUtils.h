#ifndef COMB_PERM_UTILS_H
#define COMB_PERM_UTILS_H

#include <Rcpp.h>

void SetClass(bool &IsCharacter, bool &IsLogical, 
              bool &IsInteger, bool &IsComplex, bool &IsRaw, SEXP Rv);

void SetValues(bool IsCharacter, bool IsLogical, bool &IsInteger, bool IsComplex,
               bool IsRaw, Rcpp::CharacterVector &rcppChar, std::vector<int> &vInt,
               std::vector<double> &vNum, Rcpp::ComplexVector &rcppCplx, 
               Rcpp::RawVector &rcppRaw, int &n, SEXP Rv);

void SetThreads(bool &Parallel, int maxThreads, int nRows,
                bool IsCharacter, int &nThreads, SEXP RNumThreads, int limit);

void SetRandomSample(SEXP RindexVec, SEXP RNumSamp, std::size_t &sampSize,
                     bool IsGmp, double computedRows, std::vector<double> &mySample,
                     Rcpp::Function baseSample);

std::vector<int> rleCpp(const std::vector<int> &x);
double NumPermsWithRep(const std::vector<int> &v);
double NumPermsNoRep(int n, int k);
double nChooseK(int n, int k);
double NumCombsWithRep(int n, int r);
double MultisetCombRowNumFast(int n,int r, const std::vector<int> &Reps);
double MultisetPermRowNum(int n, int r, const std::vector<int> &myReps);

// This one isn't as efficient as MultisetCombRowNumFast, however it will
// not produce negative results and is thus used in determining whether
// gmp analogs are necessary
double MultisetCombRowNum(int n, int r, const std::vector<int> &Reps);

void nextFullPerm(int *const myArray, std::size_t maxInd);
void nextPartialPerm(int *const myArray, std::size_t lastCol, std::size_t maxInd);

double GetComputedRows(bool IsMultiset, bool IsComb, bool IsRep, int n, int &m, SEXP Rm,
                       int lenFreqs, std::vector<int> &freqs, std::vector<int> &myReps);
    
#endif
