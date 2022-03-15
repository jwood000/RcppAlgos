#ifndef COMBINATIONS_REP_H
#define COMBINATIONS_REP_H

#include "NextComboSection.h"

template <typename T>
void CombinationsRep(T* mat, const std::vector<T> &v,
                     std::vector<int> &z, int n, int m, int nRows) {

    for (int count = 0, m1 = m - 1, n1 = n - 1; count < nRows;) {
        for (; z[m1] < n && count < nRows; ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                mat[count + j * nRows] = v[z[j]];
            }
        }

        nextCombSecRep(z, m1, n1);
    }
}

template <typename typeMat, typename T>
void CombinationsRep(typeMat &mat, const std::vector<T> &v,
                     std::vector<int> &z, int n,
                     int m, int strt, int nRows) {

    for (int count = strt, m1 = m - 1, n1 = n - 1; count < nRows;) {
        for (; z[m1] < n && count < nRows; ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                mat(count, j) = v[z[j]];
            }
        }

        nextCombSecRep(z, m1, n1);
    }
}

void CombinationsRep(SEXP mat, SEXP v, std::vector<int> &z,
                     int n, int m, int nRows) {

    for (int count = 0, m1 = m - 1, n1 = n - 1; count < nRows;) {
        for (; z[m1] < n && count < nRows; ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(mat, count + j * nRows, STRING_ELT(v, z[j]));
            }
        }

        nextCombSecRep(z, m1, n1);
    }
}

#endif
