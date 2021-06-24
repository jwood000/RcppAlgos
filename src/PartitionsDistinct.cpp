#include "Permutations/PermuteCount.h"
#include "Partitions/NextPartition.h"
#include "Cpp14MakeUnique.h"
#include <numeric>  // std::iota

template <typename T>
void PartsGenDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, int width, int lastElem,
                      int lastCol, int strt, int nRows) {

    int edge = 0;
    int pivot = 0;
    int tarDiff = 0;
    int boundary = 0;

    PrepareDistinctPart(z, boundary, pivot, edge, tarDiff, lastElem, lastCol);

    for (int count = strt; count < nRows; ++count,
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {

        for (int k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
    }
}

template <typename T>
void PartsGenPermDistinct(T* mat, const std::vector<T> &v,
                          std::vector<int> &z, int width,
                          int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int tarDiff = 0;
    int boundary = 0;

    PrepareDistinctPart(z, boundary, pivot, edge, tarDiff, lastElem, lastCol);

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

    for (int count = 0;;) {
        for (int j = 0, myRow = 0; j < indexRows; ++count, ++j) {
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
void PartsGenPermZeroDistinct(T* mat, const std::vector<T> &v,
                              std::vector<int> &z, int width,
                              int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int tarDiff = 0;
    int boundary = 0;

    PrepareDistinctPart(z, boundary, pivot, edge, tarDiff, lastElem, lastCol);

    for (int count = 0;;) {

        do {
            for (int k = 0; k < width; ++k) {
                mat[count + nRows * k] = v[z[k]];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
        NextDistinctGenPart(z, boundary, edge, pivot,
                            tarDiff, lastCol, lastElem);
    }
}

void PartsDistinct(int* mat, std::vector<int> &z, int width, int boundary,
                   int lastCol, int edge, int strt, int nRows) {

    for (int count = strt, tarDiff = 3; count < nRows; ++count) {
        for (int k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }

        NextDistinctPart(z, boundary, edge, tarDiff, lastCol);
    }
}

void PartsPermDistinct(int* mat, std::vector<int> &z, int width,
                       int boundary, int lastCol, int edge, int nRows) {

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

    for (int count = 0, tarDiff = 3;;
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
                           int boundary, int lastCol, int edge, int nRows) {

    int totalCount = 0;

    for (int count = 0, tarDiff = 3; z[1] == 0;
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        do {
            for (int k = 0; k < width; ++k) {
                mat[count + nRows * k] = z[k];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
        totalCount = count;
    }

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

    for (int count = totalCount, tarDiff = 3; ;
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        for (int j = 0, myRow = 0; j < indexRows; ++count, ++j) {
            for (int k = 0; k < width; ++k, ++myRow) {
                mat[count + nRows * k] = z[indexMat[myRow]];
            }
        }

        if (count >= nRows) {break;}
    }
}

template void PartsGenDistinct(int*, const std::vector<int>&,
                               std::vector<int>&, int, int, int, int, int);

template void PartsGenDistinct(double*, const std::vector<double>&,
                               std::vector<int>&, int, int, int, int, int);

template void PartsGenPermDistinct(int*, const std::vector<int>&,
                                   std::vector<int>&, int, int, int, int);

template void PartsGenPermDistinct(double*, const std::vector<double>&,
                                   std::vector<int>&, int, int, int, int);

template void PartsGenPermZeroDistinct(int*, const std::vector<int>&,
                                       std::vector<int>&, int, int, int, int);

template void PartsGenPermZeroDistinct(double*, const std::vector<double>&,
                                       std::vector<int>&, int, int, int, int);
