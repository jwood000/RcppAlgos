#pragma once

#include <vector>

using prevIterPtr = void (*const)(const std::vector<int> &freqs,
                          std::vector<int> &z,int n1, int m1);

prevIterPtr GetPrevIterPtr(bool IsComb, bool IsMult,
                           bool IsRep, bool IsFull);
