#ifndef COMB_PERM_UTILS_H
#define COMB_PERM_UTILS_H

#include <Rcpp.h>

std::vector<std::vector<int> > rleCpp(std::vector<int> &x);
double NumPermsWithRep(std::vector<int> &v);
double NumPermsNoRep(int n, int k);
double nChooseK(double n, double k);
double NumCombsWithRep(int n, int r);
double MultisetCombRowNum(int n, int r, std::vector<int> &Reps);
double MultisetPermRowNum(int n, int r, std::vector<int> &myReps);

void nextFullPerm(int *myArray, const unsigned long int &n1,
                  const unsigned long int &n2);

void nextPartialPerm(int *myArray, const unsigned long int &r,
                     const unsigned long int &r1, const unsigned long int &n,
                     const unsigned long int &lastElem);
    
#endif
