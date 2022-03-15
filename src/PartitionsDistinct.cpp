#include "Permutations/PermuteCount.h"
#include "Partitions/NextPartition.h"
#include "Cpp14MakeUnique.h"
#include "PopulateVec.h"
#include "RMatrix.h"
#include <algorithm>  // std::next_permutation
#include <numeric>    // std::iota

template <typename T>
void PartsGenDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, int width,
                      int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int initTarDf = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        initTarDf, lastElem, lastCol);

    for (int count = 0, tarDiff = initTarDf; count < nRows; ++count,
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {

        for (int k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
    }
}

template <typename T>
void PartsGenDistinct(RcppParallel::RMatrix<T> &mat,
                      const std::vector<T> &v, std::vector<int> &z,
                      int strt, int width, int lastElem,
                      int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int initTarDf = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        initTarDf, lastElem, lastCol);

    for (int count = strt, tarDiff = initTarDf; count < nRows; ++count,
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {

        for (int k = 0; k < width; ++k) {
            mat(count, k) = v[z[k]];
        }
    }
}

template <typename T>
void PartsGenPermDistinct(T* mat, const std::vector<T> &v,
                          std::vector<int> &z, int width,
                          int lastElem, int lastCol, int nRows) {
    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int initTarDf = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        initTarDf, lastElem, lastCol);

    const int indexRows = NumPermsNoRep(width, width);
    auto indexMat = FromCpp14::make_unique<int[]>(indexRows * width);

    std::vector<int> indexVec(width);
    std::iota(indexVec.begin(), indexVec.end(), 0);

    for (int i = 0, myRow = 0; i < indexRows; ++i, myRow += width) {
        for (int j = 0; j < width; ++j) {
            indexMat[myRow + j] = indexVec[j];
        }

        std::next_permutation(indexVec.begin(), indexVec.end());
    }

    for (int count = 0, tarDiff= initTarDf;;) {
        for (int j = 0, myRow = 0; j < indexRows && count < nRows; ++count, ++j) {
            for (int k = 0; k < width; ++k, ++myRow) {
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
                      std::vector<int> &z, int width,
                      int nRows, bool IsComb) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int initTarDf = 0;

    const int lastCol = width - 1;
    const int lastElem = v.size() - 1;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        initTarDf, lastElem, lastCol);

    for (int count = 0, tarDiff = initTarDf;
         edge >= 0 && (z[boundary] - z[edge]) >= tarDiff;
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {

        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
        if (count >= nRows) break;
    }

    int count = partsVec.size() / width;

    if (count < nRows) {
        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
    }
}

template <typename T>
void PartsGenPermZeroDistinct(T* mat, const std::vector<T> &v,
                              std::vector<int> &z, int width,
                              int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int initTarDf = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        initTarDf, lastElem, lastCol);

    for (int count = 0, tarDiff = initTarDf;;
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {
        do {
            for (int k = 0; k < width; ++k) {
                mat[count + nRows * k] = v[z[k]];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
    }
}

void PartsDistinct(int* mat, std::vector<int> &z, int width,
                   int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int initTarDf = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        initTarDf, lastElem, lastCol);

    for (int count = 0, tarDiff = initTarDf; count < nRows; ++count,
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        for (int k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }
    }
}

void PartsDistinct(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                   int strt, int width, int lastElem,
                   int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int initTarDf = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        initTarDf, lastElem, lastCol);

    for (int count = strt, tarDiff = initTarDf; count < nRows; ++count,
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        for (int k = 0; k < width; ++k) {
            mat(count, k) = z[k];
        }
    }
}

void PartsPermDistinct(int* mat, std::vector<int> &z, int width,
                       int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int initTarDf = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        initTarDf, lastElem, lastCol);

    const int indexRows = NumPermsNoRep(width, width);
    auto indexMat = FromCpp14::make_unique<int[]>(indexRows * width);

    std::vector<int> indexVec(width);
    std::iota(indexVec.begin(), indexVec.end(), 0);

    for (int i = 0, myRow = 0; i < indexRows; ++i, myRow += width) {
        for (int j = 0; j < width; ++j) {
            indexMat[myRow + j] = indexVec[j];
        }

        std::next_permutation(indexVec.begin(), indexVec.end());
    }

    for (int count = 0, tarDiff = initTarDf;;
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        for (int j = 0, myRow = 0; j < indexRows; ++count, ++j) {
            for (int k = 0; k < width; ++k, ++myRow) {
                mat[count + nRows * k] = z[indexMat[myRow]];
            }
        }

        if (count >= nRows) {break;}
    }
}

void PartsPermZeroDistinct(int* mat, std::vector<int> &z, int width,
                           int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;
    int initTarDf = 0;

    PrepareDistinctPart(z, boundary, pivot, edge,
                        initTarDf, lastElem, lastCol);

    for (int count = 0, tarDiff = initTarDf;;
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        do {
            for (int k = 0; k < width; ++k) {
                mat[count + nRows * k] = z[k];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
    }
}

template void PartsGenDistinct(int*, const std::vector<int>&,
                               std::vector<int>&, int, int, int, int);
template void PartsGenDistinct(double*, const std::vector<double>&,
                               std::vector<int>&, int, int, int, int);

template void PartsGenDistinct(RcppParallel::RMatrix<int> &mat,
                               const std::vector<int>&, std::vector<int>&,
                               int, int, int, int, int);
template void PartsGenDistinct(RcppParallel::RMatrix<double> &mat,
                               const std::vector<double>&, std::vector<int>&,
                               int, int, int, int, int);

template void PartsGenDistinct(std::vector<int>&, const std::vector<int>&,
                               std::vector<int>&, int, int, bool);
template void PartsGenDistinct(std::vector<double>&,
                               const std::vector<double>&,
                               std::vector<int>&, int, int, bool);

template void PartsGenPermDistinct(int*, const std::vector<int>&,
                                   std::vector<int>&, int, int, int, int);
template void PartsGenPermDistinct(double*, const std::vector<double>&,
                                   std::vector<int>&, int, int, int, int);

template void PartsGenPermZeroDistinct(int*, const std::vector<int>&,
                                       std::vector<int>&,
                                       int, int, int, int);
template void PartsGenPermZeroDistinct(double*, const std::vector<double>&,
                                       std::vector<int>&,
                                       int, int, int, int);
