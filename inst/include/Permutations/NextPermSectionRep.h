#ifndef NEXT_PERM_SECTION_REP_H
#define NEXT_PERM_SECTION_REP_H

#include <vector>

inline void NextSecRep(std::vector<int> &z, int maxInd, int lastCol) {

    for (int i = lastCol; i >= 0; --i) {
        if (z[i] != maxInd) {
            ++z[i];
            break;
        } else {
            z[i] = 0;
        }
    }
}

#endif
