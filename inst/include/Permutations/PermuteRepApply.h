#pragma once

#include "NextPermSectionRep.h"
#include "PermuteHelper.h"
#include "FunAssign.h"

template <typename T>
void PermuteRepApplyFun(SEXP res, const std::vector<T> &v,
                        SEXP vectorPass, T* ptr_vec, std::vector<int> &z,
                        int n, int m, int nRows, SEXP sexpFun, SEXP rho,
                        int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);

    for (int count = 0, maxInd = n - 1,
         lastCol = m - 1; count < nRows; ++count) {

        for (int j = 0; j < m; ++j) {
            ptr_vec[j] = v[z[j]];
        }

        FunAssign(res, vectorPass, sexpFun, rho,
                  commonType, commonLen, count, nRows, retType);
        NextSecRep(z, maxInd, lastCol);
    }
}

void PermuteRepApplyFun(SEXP res, SEXP v, SEXP vectorPass,
                        std::vector<int> &z, int n, int m,
                        int nRows, SEXP sexpFun, SEXP rho,
                        int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);

    for (int count = 0, maxInd = n - 1,
         lastCol = m - 1; count < nRows; ++count) {

        for (int j = 0; j < m; ++j) {
            SET_STRING_ELT(vectorPass, j, STRING_ELT(v, z[j]));
        }

        FunAssign(res, vectorPass, sexpFun, rho,
                  commonType, commonLen, count, nRows, retType);
        NextSecRep(z, maxInd, lastCol);
    }
}
