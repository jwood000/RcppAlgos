#pragma once

#include <vector>
#include "CppConvert/GmpxxCopy.h"

bool nextComboGroup(std::vector<int> &z, int nGrps,
                    int grpSize, int idx1, int idx2, int last1);

double numGroupCombs(int n, int numGroups, int grpSize);
void numGroupCombsGmp(mpz_class &result, int n, int numGroups, int grpSize);

std::vector<int> nthComboGroup(int n, int gSize, int r,
                               double myIndex, double total);

std::vector<int> nthComboGroupGmp(int n, int gSize, int r,
                                  const mpz_class &lowerMpz,
                                  const mpz_class &computedRowMpz);
