#ifndef PERMUTE_OPTIMIZED_H
#define PERMUTE_OPTIMIZED_H

#include <vector>

template <typename T>
void PermuteOptimized(T* mat, const std::vector<T> &v, std::vector<int> &z,
                      int n, int m, int nRows, bool IsRep);

#endif
