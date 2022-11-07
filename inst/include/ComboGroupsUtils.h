#pragma once

#include <vector>
#include <gmpxx.h>

bool nextComboGroup(std::vector<int> &z, int grpSize,
                    int idx1, int idx2, int low_one);

bool nextComboGroup(std::vector<int> &z,
                    const std::vector<int> &grpSize,
                    int idx1, int idx2, int lbound);

double numGroupCombs(int n, int numGroups, int grpSize);
void numGroupCombsGmp(mpz_class &result, int n, int numGroups, int grpSize);

std::vector<int> nthComboGroup(int n, int grpSize, int r,
                               double myIndex, double total);

std::vector<int> nthComboGroupGmp(int n, int grpSize, int r,
                                  const mpz_class &lowerMpz,
                                  const mpz_class &computedRowMpz);
