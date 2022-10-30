#pragma once

#include <vector>

using nextIterPtr = bool (*const)(const std::vector<int> &freqs,
                          std::vector<int> &z,int n1, int m1);

nextIterPtr GetNextIterPtr(bool IsComb, bool IsMult,
                           bool IsRep, bool IsFull);
