#pragma once

#include <vector>

bool keepGoing(const std::vector<int> &rpsCnt, int lastElem,
               const std::vector<int> &z, int edge, int boundary);

template <typename T>
void PartsGenMultiset(std::vector<T> &partsVec, const std::vector<T> &v,
                      const std::vector<int> &Reps, std::vector<int> &z,
                      std::size_t width, std::size_t nRows, bool IsComb);
