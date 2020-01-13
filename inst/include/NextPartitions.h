#ifndef NEXT_PARTITIONS_H
#define NEXT_PARTITIONS_H

#include <vector>

void PrepareMultiSetPart(const std::vector<int> &rpsCnt, const std::vector<int> &z,
                         int &b, int &p, int &e, int lastCol, int lastElem);

void NextDistinct(std::vector<int> &z, int &boundary,
                  int &edge, int &tarDiff, int lastCol);

void NextPartition(std::vector<int> &z, int &boundary, int &edge, int lastCol);

void NextMultiSetGenPart(std::vector<int> &rpsCnt, std::vector<int> &z,
                         int &e, int &b, int &p, int lastCol, int lastElem);

void NextRepGenPart(std::vector<int> &z, int &boundary, int &edge, 
                    int &pivot, int lastCol, int lastElem);

void NextDistinctGenPart(std::vector<int> &z, int &boundary, int &edge, 
                         int &pivot, int &tarDiff, int lastCol, int lastElem);

#endif
