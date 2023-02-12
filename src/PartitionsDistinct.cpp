#include <memory>
#include "Permutations/PermuteCount.h"
#include "Partitions/NextPartition.h"
#include "PopulateVec.h"
#include "RMatrix.h"
#include <algorithm>  // std::next_permutation
#include <numeric>    // std::iota

template <typename T>
void PartsGenDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, std::size_t width,
                      int lastElem, int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int tarDiff = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        tarDiff, lastElem, lastCol);

    for (std::size_t count = 0; count < nRows; ++count,
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
    }
}

template <typename T>
void PartsGenDistinct(RcppParallel::RMatrix<T> &mat,
                      const std::vector<T> &v, std::vector<int> &z,
                      int strt, std::size_t width, int lastElem,
                      int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int tarDiff = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        tarDiff, lastElem, lastCol);

    for (std::size_t count = strt; count < nRows; ++count,
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = v[z[k]];
        }
    }
}

template <typename T>
void PartsGenPermDistinct(T* mat, const std::vector<T> &v,
                          std::vector<int> &z, std::size_t width,
                          int lastElem, int lastCol, std::size_t nRows) {
    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int tarDiff = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        tarDiff, lastElem, lastCol);

    const int indexRows = NumPermsNoRep(width, width);
    auto indexMat = std::make_unique<int[]>(indexRows * width);

    std::vector<int> indexVec(width);
    std::iota(indexVec.begin(), indexVec.end(), 0);

    for (int i = 0, myRow = 0; i < indexRows; ++i, myRow += width) {
        for (std::size_t j = 0; j < width; ++j) {
            indexMat[myRow + j] = indexVec[j];
        }

        std::next_permutation(indexVec.begin(), indexVec.end());
    }

    for (std::size_t count = 0;;) {
        for (int j = 0, myRow = 0; j < indexRows && count < nRows;
             ++count, ++j) {
            for (std::size_t k = 0; k < width; ++k, ++myRow) {
                mat[count + nRows * k] = v[z[indexMat[myRow]]];
            }
        }

        if (count >= nRows) {break;}
        NextDistinctGenPart(z, boundary, edge, pivot,
                            tarDiff, lastCol, lastElem);
    }
}

template <typename T>
void PartsGenDistinct(std::vector<T> &partsVec, const std::vector<T> &v,
                      std::vector<int> &z, std::size_t width,
                      std::size_t nRows, bool IsComb) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int tarDiff = 0;

    const int lastCol = width - 1;
    const int lastElem = v.size() - 1;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        tarDiff, lastElem, lastCol);

    for (std::size_t count = 0; edge >= 0 &&
                                (z[boundary] - z[edge]) >= tarDiff;
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {

        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
        if (count >= nRows) break;
    }

    std::size_t count = partsVec.size() / width;

    if (count < nRows) {
        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
    }
}

template <typename T>
void PartsGenPermZeroDistinct(T* mat, const std::vector<T> &v,
                              std::vector<int> &z, std::size_t width,
                              int lastElem, int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int tarDiff = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        tarDiff, lastElem, lastCol);

    for (std::size_t count = 0;;
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {
        do {
            for (std::size_t k = 0; k < width; ++k) {
                mat[count + nRows * k] = v[z[k]];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
    }
}

void PartsDistinct(int* mat, std::vector<int> &z, std::size_t width,
                   int lastElem, int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int tarDiff = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        tarDiff, lastElem, lastCol);

    for (std::size_t count = 0; count < nRows; ++count,
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }
    }
}

void PartsDistinct(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                   int strt, std::size_t width, int lastElem,
                   int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int tarDiff = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        tarDiff, lastElem, lastCol);

    for (std::size_t count = strt; count < nRows; ++count,
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = z[k];
        }
    }
}

void PartsPermDistinct(int* mat, std::vector<int> &z, std::size_t width,
                       int lastElem, int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int tarDiff = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        tarDiff, lastElem, lastCol);

    const int indexRows = NumPermsNoRep(width, width);
    auto indexMat = std::make_unique<int[]>(indexRows * width);

    std::vector<int> indexVec(width);
    std::iota(indexVec.begin(), indexVec.end(), 0);

    for (int i = 0, myRow = 0; i < indexRows; ++i, myRow += width) {
        for (std::size_t j = 0; j < width; ++j) {
            indexMat[myRow + j] = indexVec[j];
        }

        std::next_permutation(indexVec.begin(), indexVec.end());
    }

    for (std::size_t count = 0;;
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        for (int j = 0, myRow = 0; j < indexRows; ++count, ++j) {
            for (std::size_t k = 0; k < width; ++k, ++myRow) {
                mat[count + nRows * k] = z[indexMat[myRow]];
            }
        }

        if (count >= nRows) {break;}
    }
}

void PartsPermZeroDistinct(int* mat, std::vector<int> &z, std::size_t width,
                           int lastElem, int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int tarDiff = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        tarDiff, lastElem, lastCol);

    for (std::size_t count = 0;;
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        do {
            for (std::size_t k = 0; k < width; ++k) {
                mat[count + nRows * k] = z[k];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
    }
}

template void PartsGenDistinct(int*, const std::vector<int>&,
                               std::vector<int>&, std::size_t,
                               int, int, std::size_t);
template void PartsGenDistinct(double*, const std::vector<double>&,
                               std::vector<int>&, std::size_t,
                               int, int, std::size_t);

template void PartsGenDistinct(RcppParallel::RMatrix<int> &mat,
                               const std::vector<int>&, std::vector<int>&,
                               int, std::size_t, int, int, std::size_t);
template void PartsGenDistinct(RcppParallel::RMatrix<double> &mat,
                               const std::vector<double>&, std::vector<int>&,
                               int, std::size_t, int, int, std::size_t);

template void PartsGenDistinct(std::vector<int>&, const std::vector<int>&,
                               std::vector<int>&, std::size_t,
                               std::size_t, bool);
template void PartsGenDistinct(std::vector<double>&,
                               const std::vector<double>&,
                               std::vector<int>&, std::size_t,
                               std::size_t, bool);

template void PartsGenPermDistinct(int*, const std::vector<int>&,
                                   std::vector<int>&, std::size_t,
                                   int, int, std::size_t);
template void PartsGenPermDistinct(double*, const std::vector<double>&,
                                   std::vector<int>&, std::size_t,
                                   int, int, std::size_t);

template void PartsGenPermZeroDistinct(int*, const std::vector<int>&,
                                       std::vector<int>&, std::size_t,
                                       int, int, std::size_t);
template void PartsGenPermZeroDistinct(double*, const std::vector<double>&,
                                       std::vector<int>&, std::size_t,
                                       int, int, std::size_t);
