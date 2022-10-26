#pragma once

#include "NextComboSection.h"
#include "FunAssign.h"
#include <algorithm>

template <typename T>
void MultisetComboApplyFun(SEXP res, const std::vector<T> &v,
                           SEXP vectorPass, T* ptr_vec, std::vector<int> &z,
                           int n, int m, int nRows, SEXP sexpFun, SEXP rho,
                           const std::vector<int> &freqs,
                           int commonLen = 1, int commonType = INTSXP) {

    std::vector<int> zIndex(n);
    const int retType = TYPEOF(res);

    for (int i = 0; i < n; ++i) {
        zIndex[i] = std::find(freqs.cbegin(),
                              freqs.cend(), i) - freqs.cbegin();
    }

    for (int count = 0, m1 = m - 1,
         pentExtreme = freqs.size() - m; count < nRows;) {

        for (; z[m1] < n && count < nRows; ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                ptr_vec[j] = v[z[j]];
            }

            FunAssign(res, vectorPass, sexpFun, rho,
                      commonType, commonLen, count, nRows, retType);
        }

        nextCombSecMulti(freqs, zIndex, z, m1, pentExtreme);
    }
}

void MultisetComboApplyFun(SEXP res, SEXP v, SEXP vectorPass,
                           std::vector<int> &z, int n, int m,
                           int nRows, SEXP sexpFun, SEXP rho,
                           const std::vector<int> &freqs,
                           int commonLen = 1, int commonType = INTSXP) {

    std::vector<int> zIndex(n);
    const int retType = TYPEOF(res);

    for (int i = 0; i < n; ++i) {
        zIndex[i] = std::find(freqs.cbegin(),
                              freqs.cend(), i) - freqs.cbegin();
    }

    for (int count = 0, m1 = m - 1,
         pentExtreme = freqs.size() - m; count < nRows;) {

        for (; z[m1] < n && count < nRows; ++count, ++z[m1]){
            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(vectorPass, j, STRING_ELT(v, z[j]));
            }

            FunAssign(res, vectorPass, sexpFun, rho,
                      commonType, commonLen, count, nRows, retType);
        }

        nextCombSecMulti(freqs, zIndex, z, m1, pentExtreme);
    }
}
