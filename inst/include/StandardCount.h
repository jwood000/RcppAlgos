#ifndef STANDARD_COUNT_H
#define STANDARD_COUNT_H

#include <Rcpp.h>

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
