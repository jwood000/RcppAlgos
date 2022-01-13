#ifndef POPULATE_VEC_H
#define POPULATE_VEC_H

#include <vector>

template <typename T>
void PopulateVec(const std::vector<T> &v,
                 std::vector<T> &cnstrntVec,
                 std::vector<int> &z, int &count,
                 int m, int nRows, bool IsComb);

#endif
