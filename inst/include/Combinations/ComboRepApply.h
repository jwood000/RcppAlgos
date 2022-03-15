#ifndef COMBO_REP_APPLY_H
#define COMBO_REP_APPLY_H

#include "NextComboSection.h"
#include "FunAssign.h"

template <typename T>
void ComboRepApplyFun(SEXP res, const std::vector<T> &v,
                      SEXP vectorPass, T* ptr_vec, std::vector<int> &z,
                      int n, int m, int nRows, SEXP sexpFun, SEXP rho,
                      int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);

    for (int count = 0, m1 = m - 1, n1 = n - 1; count < nRows;) {
        for (; z[m1] < n && count < nRows; ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                ptr_vec[j] = v[z[j]];
            }

            FunAssign(res, vectorPass, sexpFun, rho,
                      commonType, commonLen, count, nRows, retType);
        }

        nextCombSecRep(z, m1, n1);
    }
}

void ComboRepApplyFun(SEXP res, SEXP v, SEXP vectorPass,
                      std::vector<int> &z, int n, int m,
                      int nRows, SEXP sexpFun, SEXP rho,
                      int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);

    for (int count = 0, m1 = m - 1, n1 = n - 1; count < nRows;) {
        for (; z[m1] < n && count < nRows; ++count, ++z[m1]) {
            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(vectorPass, j, STRING_ELT(v, z[j]));
            }

            FunAssign(res, vectorPass, sexpFun, rho,
                      commonType, commonLen, count, nRows, retType);
        }

        nextCombSecRep(z, m1, n1);
    }
}

#endif
