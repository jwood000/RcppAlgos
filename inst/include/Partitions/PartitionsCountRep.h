#pragma once

double CountPartsRepLenRstrctd(
    int n, int m, const std::vector<int> &allowed, int strtLen = 0
);

double CountPartsRepLen(
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);

double CountPartsRep(
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);

double CountCompsRepLen(
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);

double CountCompsRepZNotWk(
    int n, int m, const std::vector<int> &allowed = std::vector<int>(),
    int strtLen = 0
);
