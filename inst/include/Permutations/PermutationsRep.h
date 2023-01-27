#pragma once

#include "NextPermSectionRep.h"
#include "PermuteHelper.h"
#include "RMatrix.h"

template <typename T>
void PermuteRep(T* mat, const std::vector<T> &v,
                std::vector<int> &z, std::size_t n,
                std::size_t m, std::size_t nRows) {

    for (std::size_t count = 0, maxInd = n - 1,
         lastCol = m - 1; count < nRows; ++count) {

        for (std::size_t j = 0; j < m; ++j) {
            mat[count + j * nRows] = v[z[j]];
        }

        NextSecRep(z, maxInd, lastCol);
    }
}

template <typename T>
void PermuteRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                std::vector<int> &z, std::size_t n, std::size_t m,
                std::size_t strt, std::size_t nRows) {

    for (std::size_t count = strt, maxInd = n - 1,
         lastCol = m - 1; count < nRows; ++count) {

        for (std::size_t j = 0; j < m; ++j) {
            mat(count, j) = v[z[j]];
        }

        NextSecRep(z, maxInd, lastCol);
    }
}

void PermuteRep(SEXP mat, SEXP v, std::vector<int> &z,
                std::size_t n, std::size_t m, std::size_t nRows) {

    for (std::size_t count = 0, maxInd = n - 1,
         lastCol = m - 1; count < nRows; ++count) {

        for (std::size_t j = 0; j < m; ++j) {
            SET_STRING_ELT(mat, count + j * nRows, STRING_ELT(v, z[j]));
        }

        NextSecRep(z, maxInd, lastCol);
    }
}
