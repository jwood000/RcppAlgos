#ifndef POPULATE_VEC_PERM_H
#define POPULATE_VEC_PERM_H

#include <algorithm>
#include <vector>

template <typename T>
inline void PopulateVecPerm(const std::vector<T> &v,
                            std::vector<T> &partitionsVec,
                            std::vector<int> &z, int &count,
                            int width, int maxRows) {
    do {
        for (int k = 0; k < width; ++k)
            partitionsVec.push_back(v[z[k]]);

        ++count;
    } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
}

#endif