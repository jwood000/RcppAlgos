#ifndef COMB_PERM_UTILS_H
#define COMB_PERM_UTILS_H

#include <Rcpp.h>

double NumPermsWithRep(std::vector<int> v);
double NumPermsNoRep(int n, int k);
double nChooseK(double n, double k);
double NumCombsWithRep(int n, int r);
double MultisetCombRowNum(int n, int r, std::vector<int> Reps);
double MultisetPermRowNum(int n, int r, std::vector<int> myReps);

void nextFullPerm(uint16_t *myArray, unsigned long int n1);

void nextPartialPerm(uint16_t *myArray, unsigned long int nCols, 
                     unsigned long int r1, unsigned long int r,
                     unsigned long int n1, unsigned long int n);

#endif
