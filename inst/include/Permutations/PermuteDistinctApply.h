#ifndef PERMUTE_DISTINCT_APPLY_H
#define PERMUTE_DISTINCT_APPLY_H

#include "PermuteHelper.h"
#include "FunAssign.h"

template <typename T>
void PermuteDistinctApplyFun(SEXP res, const std::vector<T> &v,
                             SEXP vectorPass, T* ptr_vec, std::vector<int> &z,
                             int n, int m, int nRows, SEXP sexpFun, SEXP rho,
                             int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);
    auto arrPerm = FromCpp14::make_unique<int[]>(n);

    for (int i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        for (int count = 0, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                ptr_vec[j] = v[arrPerm[j]];
            }

            FunAssign(res, vectorPass, sexpFun, rho,
                      commonType, commonLen, count, nRows, retType);
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                ptr_vec[j] = v[arrPerm[j]];
            }

            FunAssign(res, vectorPass, sexpFun, rho,
                      commonType, commonLen, count, nRows, retType);
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j) {
        ptr_vec[j] = v[arrPerm[j]];
    }

    FunAssign(res, vectorPass, sexpFun, rho,
              commonType, commonLen, nRows - 1, nRows, retType);
}

void PermuteDistinctApplyFun(SEXP res, SEXP v, SEXP vectorPass,
                             std::vector<int> &z, int n, int m,
                             int nRows, SEXP sexpFun, SEXP rho,
                             int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);
    auto arrPerm = FromCpp14::make_unique<int[]>(n);

    for (int i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        for (int count = 0, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(vectorPass, j, STRING_ELT(v, arrPerm[j]));
            }

            FunAssign(res, vectorPass, sexpFun, rho,
                      commonType, commonLen, count, nRows, retType);
            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(vectorPass, j, STRING_ELT(v, arrPerm[j]));
            }

            FunAssign(res, vectorPass, sexpFun, rho,
                      commonType, commonLen, count, nRows, retType);
            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j) {
        SET_STRING_ELT(vectorPass, j, STRING_ELT(v, arrPerm[j]));
    }

    FunAssign(res, vectorPass, sexpFun, rho,
              commonType, commonLen, nRows - 1, nRows, retType);
}

#endif
