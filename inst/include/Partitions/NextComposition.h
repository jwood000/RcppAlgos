#pragma once

#include <vector>

template <int one_or_zero>
void NextCompositionRep(std::vector<int> &z, int lastCol);

void NextCompositionDistinct(
    std::vector<int> &z, std::vector<int> &complement, std::vector<int> &idx,
    std::vector<int> &tailSum, int &i1, int &i2, int &myMax, int lastCol,
    int lastIdx, int target
);
