#include "Partitions/NextComposition.h"
#include "RMatrix.h"

template <typename T>
void CompsGenDistinct(
    T* mat, const std::vector<T> &v, std::vector<int> &z,
    std::vector<int> &complement, int i1, int i2, int myMax,
    int tar, std::size_t width, std::size_t nRows
) {

    // for (std::size_t count = 0, m = width - 1,
    //      q = complement.size() - 1; count < nRows; ++count,
    //      NextCompositionDistinct(z, complement, i1, i2, myMax, m, q, tar)) {
    //
    //     for (std::size_t k = 0; k < width; ++k) {
    //         mat[count + nRows * k] = v[z[k]];
    //     }
    //}
}

template <typename T>
void CompsGenDistinct(
    RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
    std::vector<int> &z, std::vector<int> &complement, int i1, int i2,
    int myMax, int tar, std::size_t strt, std::size_t width, std::size_t nRows
) {

    // for (std::size_t count = strt, m = width - 1,
    //      q = complement.size() - 1; count < nRows; ++count,
    //      NextCompositionDistinct(z, complement, i1, i2, myMax, m, q, tar)) {
    //
    //     for (std::size_t k = 0; k < width; ++k) {
    //         mat(count, k) = v[z[k]];
    //     }
    // }
}

void CompsDistinct(
    int* mat, std::vector<int> &z, std::vector<int> &complement,
    int i1, int i2, int myMax, int tar, std::size_t strt,
    std::size_t width, std::size_t nRows, std::size_t totalRows
) {

    for (std::size_t count = strt, m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow; ++count,
         NextCompositionDistinct(z, complement, i1, i2, myMax, m, q, tar)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + totalRows * k] = z[k];
        }
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat[nRows - 1 + totalRows * k] = z[k];
    }
}

void CompsDistinct(
    RcppParallel::RMatrix<int> &mat,
    std::vector<int> &z, std::vector<int> &complement, int i1, int i2,
    int myMax, int tar, std::size_t strt, std::size_t width, std::size_t nRows
) {

    // for (std::size_t count = strt, m = width - 1,
    //      q = complement.size() - 1; count < nRows; ++count,
    //     NextCompositionDistinct(z, complement, i1, i2, myMax, m, q, tar)) {
    //
    //     for (std::size_t k = 0; k < width; ++k) {
    //         mat(count, k) = z[k];
    //     }
    // }
}

// template <typename T>
// void PartsGenDistinct(std::vector<T> &partsVec, const std::vector<T> &v,
//                       std::vector<int> &z, std::size_t width,
//                       std::size_t nRows, bool IsComb) {
//
//     int edge = 0;
//     int pivot = 0;
//     int boundary = 0;
//     int tarDiff = 0;
//
//     const int lastCol = width - 1;
//     const int lastElem = v.size() - 1;
//
//     PrepareDistinctPart(z, boundary, pivot, edge,
//                         tarDiff, lastElem, lastCol);
//
//     for (std::size_t count = 0; edge >= 0 &&
//          (z[boundary] - z[edge]) >= tarDiff;
//     NextDistinctGenPart(z, boundary, edge, pivot,
//                         tarDiff, lastCol, lastElem)) {
//
//         PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
//         if (count >= nRows) break;
//     }
//
//     std::size_t count = partsVec.size() / width;
//
//     if (count < nRows) {
//         PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
//     }
// }

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
