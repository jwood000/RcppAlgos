#ifndef COMB_PERM_UTILS_H
#define COMB_PERM_UTILS_H

#include <Rcpp.h>

std::vector<int> rleCpp(const std::vector<int> &x);
double NumPermsWithRep(const std::vector<int> &v);
double NumPermsNoRep(const int n, const int k);
double nChooseK(const int n, const int k);
double NumCombsWithRep(const int n, const int r);
double MultisetCombRowNumFast(const int n, const int r, const std::vector<int> &Reps);
double MultisetPermRowNum(const int n, const int r, const std::vector<int> &myReps);

// This one isn't as efficient as MultisetCombRowNumFast, however it will
// not produce negative results and is thus used in determining whether
// gmp analogs are necessary
double MultisetCombRowNum(const int n, const int r, const std::vector<int> &Reps);
    
void nextFullPerm(int *myArray, const unsigned long int n1,
                  const unsigned long int n2);

void nextPartialPerm(int *myArray, const unsigned long int r,
                     const unsigned long int r1, const unsigned long int n,
                     const unsigned long int lastElem);
    
#endif
