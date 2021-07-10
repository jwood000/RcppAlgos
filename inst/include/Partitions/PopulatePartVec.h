#ifndef POPULATE_PART_VEC_H
#define POPULATE_PART_VEC_H

#include <algorithm>
#include <vector>

template <typename T>
inline void PopulatePartVec(const std::vector<T> &v,
                            std::vector<T> &partsVec,
                            std::vector<int> &z, int &count,
                            int width, int nRows, bool IsComb) {

    if (IsComb) {
        for (int k = 0; k < width; ++k) {
            partsVec.push_back(v[z[k]]);
        }

        ++count;
    } else {
        do {
            for (int k = 0; k < width; ++k) {
                partsVec.push_back(v[z[k]]);
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);
    }
}

#endif