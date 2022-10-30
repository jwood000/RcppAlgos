#pragma once

#include "Partitions/PartitionsTypes.h"

void GetLastPart(int* mat, std::vector<int> &z, int m, int nRows);

void PrepareMultisetPart(std::vector<int> &rpsCnt,
                         const std::vector<int> &z, int &b,
                         int &p, int &e, int lastCol, int lastElem);

void PrepareRepPart(const std::vector<int> &z, int &boundary,
                    int &pivot, int &edge, int lastElem, int lastCol);

void PrepareDistinctPart(const std::vector<int> &z, int &boundary,
                         int &pivot, int &edge, int &tarDiff,
                         int lastElem, int lastCol);

void NextDistinctPart(std::vector<int> &z, int &boundary,
                      int &edge, int &tarDiff, int lastCol);

void NextRepPart(std::vector<int> &z, int &boundary, int &edge, int lastCol);

void NextMultisetGenPart(std::vector<int> &rpsCnt, std::vector<int> &z,
                         int &e, int &b, int &p, int lastCol, int lastElem);

void NextRepGenPart(std::vector<int> &z, int &boundary, int &edge,
                    int &pivot, int lastCol, int lastElem);

void NextDistinctGenPart(std::vector<int> &z, int &boundary, int &edge,
                         int &pivot, int &tarDiff, int lastCol, int lastElem);

using nextPartsPtr = void (*const)(std::vector<int> &rpsCnt,
                           std::vector<int> &z, int &e, int &b, int &p,
                           int &tarDiff, int lastCol, int lastElem);

nextPartsPtr GetNextPartsPtr(PartitionType ptype, bool IsGen,
                             bool IsComp);
