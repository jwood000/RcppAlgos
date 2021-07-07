#ifndef NEXT_COMBO_SECTION_H
#define NEXT_COMBO_SECTION_H

#include <vector>

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

inline void nextCombSec(std::vector<int> &z, int m1, int nMinusM) {

    for (int i = m1 - 1; i >= 0; --i) {
        if (z[i] != (nMinusM + i)) {
            ++z[i];

            for (int j = i; j < m1; ++j)
                z[j + 1] = z[j] + 1;

            break;
        }
    }
}

inline void nextCombSecRep(std::vector<int> &z, int m1, int n1) {

    for (int i = m1 - 1; i >= 0; --i) {
        if (z[i] != n1) {
            ++z[i];

            for (int j = i, val = z[i]; j < m1; ++j)
                z[j + 1] = val;

            break;
        }
    }
}

inline void nextCombSecMulti(const std::vector<int> &freqs,
                             const std::vector<int> &zIndex,
                             std::vector<int> &z, int m1,
                             int pentExtreme) {

    for (int i = m1 - 1; i >= 0; --i) {
        if (z[i] != freqs[pentExtreme + i]) {
            ++z[i];

            for (int j = i + 1, k = zIndex[z[i]] + 1; j <= m1; ++j, ++k)
                z[j] = freqs[k];

            break;
        }
    }
}

#endif
