#ifndef COMBINATIONS_MULTISET_H
#define COMBINATIONS_MULTISET_H

#include "NextComboSection.h"
#include <algorithm>  // std::find

template <typename T>
void MultisetCombination(T* mat, const std::vector<T> &v,
                         std::vector<int> &z, int n, int m,
                         int nRows, const std::vector<int> &freqs) {

    std::vector<int> zIndex(n);

    for (int i = 0; i < n; ++i) {
        zIndex[i] = std::find(freqs.cbegin(),
                              freqs.cend(), i) - freqs.cbegin();
    }

    // pentExtreme is the location in freqs that represents
    // the maximal value of the second to the last element

    for (int count = 0, m1 = m - 1,
         pentExtreme = freqs.size() - m; count < nRows;) {

        for (; z[m1] < n && count < nRows; ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                mat[count + j * nRows] = v[z[j]];
            }
        }

        nextCombSecMulti(freqs, zIndex, z, m1, pentExtreme);
    }
}

template <typename typeMat, typename T>
void MultisetCombination(typeMat &mat, const std::vector<T> &v,
                         std::vector<int> &z, int n, int m, int strt,
                         int nRows, const std::vector<int> &freqs) {

    std::vector<int> zIndex(n);

    for (int i = 0; i < n; ++i) {
        zIndex[i] = std::find(freqs.cbegin(),
                              freqs.cend(), i) - freqs.cbegin();
    }

    // pentExtreme is the location in freqs that represents
    // the maximal value of the second to the last element

    for (int count = strt, m1 = m - 1,
         pentExtreme = freqs.size() - m; count < nRows;) {

        for (; z[m1] < n && count < nRows; ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                mat(count, j) = v[z[j]];
            }
        }

        nextCombSecMulti(freqs, zIndex, z, m1, pentExtreme);
    }
}

void MultisetCombination(SEXP mat, SEXP v, std::vector<int> &z,
                         int n, int m, int nRows,
                         const std::vector<int> &freqs) {

    std::vector<int> zIndex(n);

    for (int i = 0; i < n; ++i) {
        zIndex[i] = std::find(freqs.cbegin(),
                              freqs.cend(), i) - freqs.cbegin();
    }

    // pentExtreme is the location in freqs that represents
    // the maximal value of the second to the last element

    for (int count = 0, m1 = m - 1,
         pentExtreme = freqs.size() - m; count < nRows;) {

        for (; z[m1] < n && count < nRows; ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(mat, count + j * nRows, STRING_ELT(v, z[j]));
            }
        }

        nextCombSecMulti(freqs, zIndex, z, m1, pentExtreme);
    }
}

#endif
