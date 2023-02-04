#pragma once

#include "PermuteHelper.h"
#include "RMatrix.h"

template <typename T>
void PermuteMultiset(T* mat, const std::vector<T> &v,
                     std::vector<int> &z, std::size_t n, std::size_t m,
                     std::size_t nRows, const std::vector<int> &freqs) {

    const std::size_t lenFreqs = z.size();
    auto arrPerm = std::make_unique<int[]>(lenFreqs);

    for (std::size_t i = 0; i < lenFreqs; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == lenFreqs) {
        for (std::size_t count = 0, numR1 = nRows - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (std::size_t j = 0; j < m; ++j) {
                mat[count + j * nRows] = v[arrPerm[j]];
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (std::size_t count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

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
void PermuteMultiset(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     std::vector<int> &z, std::size_t n, std::size_t m,
                     std::size_t strt, std::size_t nRows,
                     const std::vector<int> &freqs) {

    const std::size_t lenFreqs = z.size();
    auto arrPerm = std::make_unique<int[]>(lenFreqs);

    for (std::size_t i = 0; i < lenFreqs; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == lenFreqs) {
        for (std::size_t count = strt, numR1 = nRows - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (std::size_t j = 0; j < m; ++j) {
                mat(count, j) = v[arrPerm[j]];
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (std::size_t count = strt, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

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

void PermuteMultiset(SEXP mat, SEXP v, std::vector<int> &z,
                     std::size_t n, std::size_t m, std::size_t nRows,
                     const std::vector<int> &freqs) {

    const std::size_t lenFreqs = z.size();
    auto arrPerm = std::make_unique<int[]>(lenFreqs);

    for (std::size_t i = 0; i < lenFreqs; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == lenFreqs) {
        for (std::size_t count = 0, numR1 = nRows - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (std::size_t j = 0; j < m; ++j) {
                SET_STRING_ELT(mat, count + j * nRows,
                               STRING_ELT(v, arrPerm[j]));
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (std::size_t count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

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
