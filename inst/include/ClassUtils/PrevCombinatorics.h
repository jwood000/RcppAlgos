#ifndef PREV_COMBINATORICS_H
#define PREV_COMBINATORICS_H

#include <vector>

using prevIterPtr = void (*const)(const std::vector<int> &freqs,
                          std::vector<int> &z,int n1, int m1);

prevIterPtr GetPrevIterPtr(bool IsComb, bool IsMult,
                           bool IsRep, bool IsFull);

#endif
