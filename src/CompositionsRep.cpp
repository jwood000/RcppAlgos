#include "Partitions/NextComposition.h"
#include "RMatrix.h"

template <int one_or_zero, typename T>
int CompsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                std::size_t width, std::size_t nRows) {

    for (std::size_t count = 0, lastCol = width - 1; count < nRows; ++count,
         NextCompositionRep<one_or_zero>(z, lastCol)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
    }

    return 1;
}

template <int one_or_zero, typename T>
int CompsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                std::vector<int> &z, std::size_t strt,
                std::size_t width, std::size_t nRows) {

    for (std::size_t count = strt, lastCol = width - 1; count < nRows; ++count,
         NextCompositionRep<one_or_zero>(z, lastCol)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = v[z[k]];
        }
    }

    return 1;
}

template <int one_or_zero>
int CompsRep(int* mat, std::vector<int> &z,
             std::size_t width, std::size_t nRows) {

    for (std::size_t count = 0, lastCol = width - 1; count < nRows; ++count,
         NextCompositionRep<one_or_zero>(z, lastCol)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }
    }

    return 1;
}

template <int one_or_zero>
int CompsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
             std::size_t strt, std::size_t width, std::size_t nRows) {

    for (std::size_t count = strt, lastCol = width - 1; count < nRows; ++count,
         NextCompositionRep<one_or_zero>(z, lastCol)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = z[k];
        }
    }

    return 1;
}

template int CompsGenRep<0>(int*, const std::vector<int>&,
                            std::vector<int>&, std::size_t, std::size_t);
template int CompsGenRep<0>(double*, const std::vector<double>&,
                            std::vector<int>&, std::size_t, std::size_t);
template int CompsGenRep<1>(int*, const std::vector<int>&,
                            std::vector<int>&, std::size_t, std::size_t);
template int CompsGenRep<1>(double*, const std::vector<double>&,
                            std::vector<int>&, std::size_t, std::size_t);

template int CompsGenRep<0>(RcppParallel::RMatrix<int>&,
                            const std::vector<int>&, std::vector<int>&,
                            std::size_t, std::size_t, std::size_t);
template int CompsGenRep<0>(RcppParallel::RMatrix<double>&,
                            const std::vector<double>&, std::vector<int>&,
                            std::size_t, std::size_t, std::size_t);
template int CompsGenRep<1>(RcppParallel::RMatrix<int>&,
                            const std::vector<int>&, std::vector<int>&,
                            std::size_t, std::size_t, std::size_t);
template int CompsGenRep<1>(RcppParallel::RMatrix<double>&,
                            const std::vector<double>&, std::vector<int>&,
                            std::size_t, std::size_t, std::size_t);

template int CompsRep<0>(int* mat, std::vector<int>&,
                         std::size_t, std::size_t);
template int CompsRep<1>(int* mat, std::vector<int>&,
                         std::size_t, std::size_t);

template int CompsRep<0>(RcppParallel::RMatrix<int>&, std::vector<int>&,
                         std::size_t, std::size_t, std::size_t);
template int CompsRep<1>(RcppParallel::RMatrix<int>&, std::vector<int>&,
                         std::size_t, std::size_t, std::size_t);
