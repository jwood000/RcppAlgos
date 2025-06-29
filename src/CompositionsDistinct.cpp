#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/NextComposition.h"
#include "RMatrix.h"
#include <algorithm> // std::set_difference
#include <numeric>

std::vector<int> PrepareComplement(const std::vector<int> &z) {

    std::vector<int> complement;
    std::vector<int> myRange(z.back());
    std::iota(myRange.begin(), myRange.end(), 1);
    std::set_difference(myRange.begin(), myRange.end(), z.begin(), z.end(),
                        std::inserter(complement, complement.begin()));
    return complement;
}

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

void CompsDistinctWorker(
    int* mat, std::vector<int> &z, std::vector<int> &complement,
    int i1, int i2, int myMax, int tar, std::size_t strt,
    std::size_t width, std::size_t nRows, std::size_t totalRows
) {

    // If the size of complement is 0 then we know we are simply taking
    // permutations of z. Same thing with size 1 but a little reasoning is
    // needed. In NextCompositionDistinct, we are swapping the last two
    // elements in z with elements in complement. If complement only has one
    // element, we are reduced to the case of taking permutations of z.
    if (complement.size() < 2) {
        std::size_t count = strt;

        do {
            for (std::size_t k = 0; k < width; ++k) {
                mat[count + totalRows * k] = z[k];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        return;
    }

    std::vector<int> idx;
    std::vector<int> tailSum;

    for (std::size_t count = strt, m = width - 1, q = complement.size() - 1,
         lastRow = nRows - 1; count < lastRow; ++count) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + totalRows * k] = z[k];
        }

        NextCompositionDistinct(
            z, complement, idx, tailSum, i1, i2, myMax, m, q, tar
        );
    }

    for (std::size_t k = 0; k < width; ++k) {
        mat[nRows - 1 + totalRows * k] = z[k];
    }
}

void CompsDistinct(int* mat, std::vector<int> &z,
                   std::size_t width, std::size_t nRows) {

    std::vector<int> complement = PrepareComplement(z);

    int strt = 0;
    const int tar = std::accumulate(z.cbegin(), z.cend(), 0);
    const int nz  = std::count(z.cbegin(), z.cend(), 0);

    for (int i = width - nz, j = nz, nextStep = 0; i < width; ++i, --j) {

        nextStep += CountCompsDistinctLen(tar, i);

        CompsDistinctWorker(
            mat, z, complement, 0, complement.size() - 1,
            z.back(), tar, strt, width, nextStep, nRows
        );

        strt = nextStep;

        std::iota(z.begin() + j - 1, z.end(), 1);
        z.back() = tar - static_cast<int>((i * (i + 1)) / 2);
        complement = PrepareComplement(z);
    }

    CompsDistinctWorker(
        mat, z, complement, 0, complement.size() - 1,
        z.back(), tar, strt, width, nRows, nRows
    );
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
