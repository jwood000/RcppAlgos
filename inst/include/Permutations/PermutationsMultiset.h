#ifndef PERMUTATIONS_MULTISET_H
#define PERMUTATIONS_MULTISET_H

#include "PermuteHelper.h"
#include "RMatrix.h"

template <typename T>
void PermuteMultiset(T* mat, const std::vector<T> &v,
                     std::vector<int> &z, int n, int m,
                     int nRows, const std::vector<int> &freqs) {

    const int lenFreqs = z.size();
    auto arrPerm = FromCpp14::make_unique<int[]>(lenFreqs);

    for (int i = 0; i < lenFreqs; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == lenFreqs) {
        for (int count = 0, numR1 = nRows - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat[count + j * nRows] = v[arrPerm[j]];
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

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
void PermuteMultiset(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                     std::vector<int> &z, int n, int m, int strt,
                     int nRows, const std::vector<int> &freqs) {

    const int lenFreqs = z.size();
    auto arrPerm = FromCpp14::make_unique<int[]>(lenFreqs);

    for (int i = 0; i < lenFreqs; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == lenFreqs) {
        for (int count = strt, numR1 = nRows - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                mat(count, j) = v[arrPerm[j]];
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = strt, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

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

void PermuteMultiset(SEXP mat, SEXP v, std::vector<int> &z, int n,
                     int m, int nRows, const std::vector<int> &freqs) {

    const int lenFreqs = z.size();
    auto arrPerm = FromCpp14::make_unique<int[]>(lenFreqs);

    for (int i = 0; i < lenFreqs; ++i) {
        arrPerm[i] = z[i];
    }

    if (m == lenFreqs) {
        for (int count = 0, numR1 = nRows - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(mat, count + j * nRows,
                               STRING_ELT(v, arrPerm[j]));
            }

            nextFullPerm(arrPerm.get(), maxInd);
        }
    } else {
        for (int count = 0, numR1 = nRows - 1, lastCol = m - 1,
             maxInd = lenFreqs - 1; count < numR1; ++count) {

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
