#pragma once

#include <vector>

bool IsComplementZeroBased(bool firstZero, bool isWeak, bool IsGen);

template <int one_or_zero>
void NextCompositionRep(std::vector<int> &z, int lastCol);

void NextCompositionDistinct(
    std::vector<int> &z, std::vector<int> &complement, std::vector<int> &idx,
    std::vector<int> &tailSum, int &i1, int &i2, int &myMax, int lastCol,
    int lastIdx, int target
);

int CompsDistinctSetup(
    const std::vector<int> &z, std::vector<int> &complement,
    int &tar, int &idx_1, int &idx_2, int &myMax, int cap, bool startAtZero
);
