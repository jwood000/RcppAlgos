#pragma once

#include <vector>

void UpdateAllowed(
    std::vector<char> &mask, std::vector<int> &allowed, int i,
    int new_val, int width, int n, int cur_val, int partial_sum
);

double CountPartsDistLenRstrctd(
    int n, int m, const std::vector<int> &allowed, int strtLen = 0
);

double CountPartsDistinctLen(
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);

double CountPartsDistinct(
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);

double CountPartsDistinctMultiZero(
    int n, int m, const std::vector<int> &allowed, int strtLen
);

double CountPartsDistinctRstrctdMZ(
    int n, int m, const std::vector<int> &allowed, int strtLen
);

double CountCompDistLenRstrctd(
    int n, int m, const std::vector<int> &allowed, int strtLen = 0
);

double CountPartsPermDistinctRstrctdMZ(
    int n, int m, const std::vector<int> &allowed, int strtLen
);

double CountCompsDistinctLen(
    int n, int m, const std::vector<int> &allowed = std::vector<int>(), int strtLen = 0
);

double CountCompsDistinctMultiZero(
    int n, int m, const std::vector<int> &allowed, int strtLen
);

double CountCompsDistinctMZWeak(
    int n, int m, const std::vector<int> &allowed, int strtLen
);
