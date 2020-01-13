#ifndef COMBO_GROUPS_UTILS_H
#define COMBO_GROUPS_UTILS_H

#include "GmpDependUtils.h" // Needed in ComboGroups.cpp.

bool nextComboGroup(std::vector<int> &z, int nGrps, 
                    int grpSize, int idx1, int idx2, int last1);

double numGroupCombs(int n, int numGroups, int grpSize);
void numGroupCombsGmp(mpz_t result, int n, int numGroups, int grpSize);

std::vector<int> nthComboGroup(int n, int gSize, int r,
                               double myIndex, double total);
    
std::vector<int> nthComboGroupGmp(int n, int gSize, int r,
                                  mpz_t lowerMpz, mpz_t computedRowMpz);

#endif
