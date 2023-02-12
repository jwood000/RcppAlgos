#pragma once

#include "PermuteHelper.h"
#include "RMatrix.h"

template <typename T>
void PermuteDistinct(T* mat, const std::vector<T> &v,
                     std::vector<int> &z, std::size_t n,
                     std::size_t m, std::size_t nRows) {

    auto arrPerm = std::make_unique<int[]>(n);

    for (std::size_t i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        for (std::size_t count = 0, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (std::size_t j = 0; j < m; ++j) {
                mat[count + j * nRows] = v[arrPerm[j]];
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (std::size_t count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (std::size_t j = 0; j < m; ++j) {
                mat[count + j * nRows] = v[arrPerm[j]];
            }

            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (std::size_t j = 0; j < m; ++j) {
        mat[nRows - 1 + j * nRows] = v[arrPerm[j]];
    }
}

template <typename T>
void PermuteDistinct(RcppParallel::RMatrix<T> &mat,
                     const std::vector<T> &v, std::vector<int> &z,
                     std::size_t n, std::size_t m,
                     std::size_t strt, std::size_t nRows) {

    auto arrPerm = std::make_unique<int[]>(n);

    for (std::size_t i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        for (std::size_t count = strt, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (std::size_t j = 0; j < m; ++j) {
                mat(count, j) = v[arrPerm[j]];
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (std::size_t count = strt, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (std::size_t j = 0; j < m; ++j) {
                mat(count, j) = v[arrPerm[j]];
            }

            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (std::size_t j = 0; j < m; ++j) {
        mat(nRows - 1, j) = v[arrPerm[j]];
    }
}

void PermuteDistinct(SEXP mat, SEXP v,
                     std::vector<int> &z, std::size_t n,
                     std::size_t m, std::size_t nRows) {

    auto arrPerm = std::make_unique<int[]>(n);

    for (std::size_t i = 0; i < n; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == n) {
        for (std::size_t count = 0, numR1 = nRows - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (std::size_t j = 0; j < m; ++j) {
                SET_STRING_ELT(mat, count + j * nRows,
                               STRING_ELT(v, arrPerm[j]));
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (std::size_t count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = n - 1; count < numR1; ++count) {

            for (std::size_t j = 0; j < m; ++j) {
                SET_STRING_ELT(mat, count + j * nRows,
                               STRING_ELT(v, arrPerm[j]));
            }

            nextPartialPerm(arrPerm.get(), lastCol, maxInd);
        }
    }

    // Get last permutation
    for (std::size_t j = 0; j < m; ++j) {
        SET_STRING_ELT(mat, nRows - 1 + j * nRows,
                       STRING_ELT(v, arrPerm[j]));
    }
}
