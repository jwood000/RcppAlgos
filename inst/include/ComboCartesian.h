#pragma once

#include <vector>

void comboGrid(std::vector<int> &cartCombs,
               std::vector<int> &lastCol,
               std::vector<int> &lenGrps,
               const std::vector<std::vector<int>> &myVecs,
               const std::vector<int> &primes, bool IsRep);
