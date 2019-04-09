#ifndef COMB_PERM_UTILS_H
#define COMB_PERM_UTILS_H

#include <Rcpp.h>

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
    
void nextFullPerm(int *myArray, unsigned long int n1,
                  unsigned long int n2);

void nextPartialPerm(int *myArray, unsigned long int r,
                     unsigned long int r1, unsigned long int n,
                     unsigned long int lastElem);
    
#endif
