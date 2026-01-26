#pragma once

#include <vector>

bool IsMaximizedGreedySuffix(const std::vector<int>& v, int target, int nz);

bool IsComplementZeroBased(bool firstZero, bool isWeak, bool IsGen);

int GetFirstPartitionDistinct(const std::vector<int> &v, std::vector<int> &z,
                              int target, int m, int lenV);

int NextDistinctBlock(const std::vector<int> &v, std::vector<int> &idx,
                      std::vector<int> &tailSum, int target, int m);

int NextDistinctBlock2(const std::vector<int> &v, std::vector<int> &idx,
                       int target, int maxLast);

int CompsDistinctSetup(
    const std::vector<int> &z, std::vector<int> &complement,
    int &tar, int &idx_1, int &idx_2, int &myMax, int idx_max,
    bool startAtZero, int zeroBudget
);
