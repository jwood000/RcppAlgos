#pragma once

#include <vector>

template <typename T>
void PermuteOptimized(T* mat, const std::vector<T> &v, std::vector<int> &z,
                      std::size_t n, std::size_t m, std::size_t nRows,
                      bool IsRep);
