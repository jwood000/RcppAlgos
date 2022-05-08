#ifndef STD_COMBO_COUNT_H
#define STD_COMBO_COUNT_H

#include <vector>

double nChooseK(int n, int m);
double NumCombsWithRep(int n, int m);
double MultisetCombRowNumFast(int n,int m, const std::vector<int> &Reps);

// This one isn't as efficient as MultisetCombRowNumFast, however it will
// not produce negative results and is thus used in determining whether
// gmp analogs are necessary
double MultisetCombRowNum(int n, int m, const std::vector<int> &Reps);

#endif
