#pragma once

#include <vector>

template <typename T>
void PopulateVec(const std::vector<T> &v,
                 std::vector<T> &cnstrntVec,
                 std::vector<int> &z, std::size_t &count,
                 std::size_t m, std::size_t nRows, bool IsComb);
