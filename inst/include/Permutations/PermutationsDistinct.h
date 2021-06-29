#ifndef PERMUTATIONS_DISTINCT_H
#define PERMUTATIONS_DISTINCT_H

#include "PermuteHelper.h"
#include "RMatrix.h"

template <typename T>
void PermuteDistinct(T* mat, const std::vector<T> &v,
                     std::vector<int> &z, int n, int m, int nRows) {

    auto arrPerm = FromCpp14::make_unique<int[]>(n);

    for (int i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        for (int count = 0, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat[count + j * nRows] = v[arrPerm[j]];
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat[count + j * nRows] = v[arrPerm[j]];
            }

            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j) {
        mat[nRows - 1 + j * nRows] = v[arrPerm[j]];
    }
}

template <typename T>
void PermuteDistinct(RcppParallel::RMatrix<T> &mat,
                     const std::vector<T> &v, std::vector<int> &z,
                     int n, int m, int strt, int nRows) {

    auto arrPerm = FromCpp14::make_unique<int[]>(n);

    for (int i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        for (int count = strt, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat(count, j) = v[arrPerm[j]];
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = strt, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat(count, j) = v[arrPerm[j]];
            }

            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j) {
        mat(nRows - 1, j) = v[arrPerm[j]];
    }
}

void PermuteDistinct(SEXP mat, SEXP v,
                     std::vector<int> &z, int n,
                     int m, int nRows) {

    auto arrPerm = FromCpp14::make_unique<int[]>(n);

    for (int i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        for (int count = 0, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(mat, count + j * nRows,
                               STRING_ELT(v, arrPerm[j]));
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(mat, count + j * nRows,
                               STRING_ELT(v, arrPerm[j]));
            }

            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (int j = 0; j < m; ++j) {
        SET_STRING_ELT(mat, nRows - 1 + j * nRows,
                       STRING_ELT(v, arrPerm[j]));
    }
}

#endif
