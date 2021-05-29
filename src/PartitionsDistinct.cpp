#include "Partitions/PopulateVecPerm.h"
#include "Permutations/PermuteCount.h"
#include "Partitions/NextPartition.h"
#include "Cpp14MakeUnique.h"
#include <numeric>  // std::iota

template <typename T>
void PartsGenDistinct(T* mat, const std::vector<T> &v,
                      std::vector<int> &z, int m, int lastElem,
                      int lastCol, int strt, int nRows) {

    int edge = 0;
    int pivot = 0;
    int tarDiff = 0;
    int boundary = 0;
    PrepareDistinctPart(z, boundary, pivot, edge, tarDiff, lastElem, lastCol);

    for (int count = strt; count < nRows; ++count,
         NextDistinctGenPart(z, boundary, edge, pivot,
                             tarDiff, lastCol, lastElem)) {

        for (int k = 0; k < m; ++k)
            mat[count + nRows * k] = v[z[k]];
    }
}

template <typename T>
void PartsGenPermDistinct(std::vector<T> &partitionsVec,
                          const std::vector<T> &v, std::vector<int> &z,
                          int m, int lastElem, int lastCol, int maxRows) {

    int edge = 0;
    int pivot = 0;
    int count = 0;
    int tarDiff = 0;
    int boundary = 0;
    PrepareDistinctPart(z, boundary, pivot, edge, tarDiff, lastElem, lastCol);

    while (edge >= 0 && (z[boundary] - z[edge]) >= tarDiff) {
        PopulateVecPerm(v, partitionsVec, z, count, m, maxRows);

        if (count >= maxRows)
            break;

        NextDistinctGenPart(z, boundary, edge, pivot,
                            tarDiff, lastCol, lastElem);
    }

    if (count < maxRows)
        PopulateVecPerm(v, partitionsVec, z, count, m, maxRows);
}

void PartsDistinct(int* mat, std::vector<int> &z,
                   int m, int boundary, int lastCol,
                   int edge, int strt, int nRows) {

    int tarDiff = 3;

    for (int count = strt; count < nRows; ++count) {
        for (int k = lastCol; k >= 0 && z[k]; --k)
            mat[count + nRows * k] = z[k];

        NextDistinctPart(z, boundary, edge, tarDiff, lastCol);
    }
}

// mIsNull && IncludeZero
void PartsPermDistinct(int* mat, std::vector<int> &z,
                       int m, int boundary, int lastCol,
                       int edge, int nRows) {
    int tarDiff = 3;

    for (int count = 0; ;
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        // index of first non-zero element
        const auto it = std::find_if(z.cbegin(), z.cend(),
                                     [](int z_i) {return z_i;});

        const int nz = std::distance(z.cbegin(), it);

        do {
            for (int k = lastCol; k >= 0 && z[k]; --k)
                mat[count + nRows * k] = z[k];

            ++count;
        } while (std::next_permutation(z.begin() + nz, z.end()) &&
                 count < nRows);

        if (count >= nRows) {break;}
    }
}

// !mIsNull || !IncludeZero
void PartsLenPermDistinct(int* mat, std::vector<int> &z,
                          int m, int boundary, int lastCol,
                          int edge, int nRows) {
    int tarDiff = 3;

    const int indexRows = NumPermsNoRep(m, m);
    auto indexMat = FromCpp14::make_unique<int[]>(indexRows * m);

    std::vector<int> indexVec(m);
    std::iota(indexVec.begin(), indexVec.end(), 0);

    for (int i = 0, myRow = 0; i < indexRows; ++i, myRow += m) {
        for (int j = 0; j < m; ++j)
            indexMat[myRow + j] = indexVec[j];

        std::next_permutation(indexVec.begin(), indexVec.end());
    }

    for (int count = 0; ;
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        for (int j = 0, myRow = 0; j < indexRows; ++count, ++j)
            for (int k = 0; k < m; ++k, ++myRow)
                mat[count + nRows * k] = z[indexMat[myRow]];

        if (count >= nRows) {break;}
    }
}

template void PartsGenDistinct(int*, const std::vector<int>&,
                               std::vector<int>&, int, int, int, int, int);

template void PartsGenDistinct(double*, const std::vector<double>&,
                               std::vector<int>&, int, int, int, int, int);

template void PartsGenPermDistinct(std::vector<int>&,
                                   const std::vector<int>&,
                                   std::vector<int>&, int, int, int, int);

template void PartsGenPermDistinct(std::vector<double>&,
                                   const std::vector<double>&,
                                   std::vector<int>&, int, int, int, int);
