#pragma once

#include <vector>

double CountPartsMultiset(
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen = 0
);

double CountCompsMultiset(
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen = 0
);

double CountCompsMultisetZero(
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
);

double CountCompsMultisetWeak(
    int n, int m, const std::vector<int>& allowed,
    const std::vector<int>& Reps, int strtLen
);
