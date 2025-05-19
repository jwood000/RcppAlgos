#include "Partitions/NextComposition.h"
#include "RMatrix.h"

template <typename T>
void CompsGenDistinct(
    T* mat, const std::vector<T> &v, std::vector<int> &z,
    std::vector<int> &complement, int i1, int i2, int myMax,
    int tar, std::size_t width, std::size_t nRows
) {

    for (std::size_t count = 0, m = width - 1,
         q = complement.size() - 1; count < nRows; ++count,
         NextCompositionDistinct(z, complement, i1, i2, myMax, m, q, tar)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
    }
}

template <typename T>
void CompsGenDistinct(
    RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
    std::vector<int> &z, std::vector<int> &complement, int i1, int i2,
    int myMax, int tar, std::size_t strt, std::size_t width, std::size_t nRows
) {

    for (std::size_t count = strt, m = width - 1,
         q = complement.size() - 1; count < nRows; ++count,
         NextCompositionDistinct(z, complement, i1, i2, myMax, m, q, tar)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = v[z[k]];
        }
    }
}

void CompsDistinct(
    int* mat, std::vector<int> &z, std::vector<int> &complement, int i1,
    int i2, int myMax, int tar, std::size_t width, std::size_t nRows
) {

    for (std::size_t count = 0, m = width - 1,
         q = complement.size() - 1; count < nRows; ++count,
         NextCompositionDistinct(z, complement, i1, i2, myMax, m, q, tar)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }
    }
}

void CompsDistinct(
    RcppParallel::RMatrix<int> &mat,
    std::vector<int> &z, std::vector<int> &complement, int i1, int i2,
    int myMax, int tar, std::size_t strt, std::size_t width, std::size_t nRows
) {

    for (std::size_t count = strt, m = width - 1,
         q = complement.size() - 1; count < nRows; ++count,
        NextCompositionDistinct(z, complement, i1, i2, myMax, m, q, tar)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = z[k];
        }
    }
}

template void CompsGenDistinct(int*, const std::vector<int>&,
                               std::vector<int>&, std::vector<int>&, int,
                               int, int, int, std::size_t, std::size_t);
template void CompsGenDistinct(double*, const std::vector<double>&,
                               std::vector<int>&, std::vector<int>&, int,
                               int, int, int, std::size_t, std::size_t);

template void CompsGenDistinct(
    RcppParallel::RMatrix<int>&, const std::vector<int>&, std::vector<int>&,
    std::vector<int>&, int, int, int, int, std::size_t, std::size_t, std::size_t
);

template void CompsGenDistinct(
    RcppParallel::RMatrix<double>&, const std::vector<double>&,
    std::vector<int>&, std::vector<int>&, int, int, int, int,
    std::size_t, std::size_t, std::size_t
);
