#include "Partitions/PopulateVecPerm.h"
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
void PartsGenPermDistinct(std::vector<T> &partitionsVec,
                          const std::vector<T> &v, std::vector<int> &z,
                          int width, int lastElem, int lastCol, int maxRows) {

    int edge = 0;
    int pivot = 0;
    int count = 0;
    int tarDiff = 0;
    int boundary = 0;
    
    PrepareDistinctPart(z, boundary, pivot, edge, tarDiff, lastElem, lastCol);

    while (edge >= 0 && (z[boundary] - z[edge]) >= tarDiff) {
        PopulateVecPerm(v, partitionsVec, z, count, width, maxRows);

        if (count >= maxRows)
            break;

        NextDistinctGenPart(z, boundary, edge, pivot,
                            tarDiff, lastCol, lastElem);
    }

    if (count < maxRows) {
        PopulateVecPerm(v, partitionsVec, z, count, width, maxRows);
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

// mIsNull && IncludeZero
void PartsPermDistinct(int* mat, std::vector<int> &z, int width,
                       int boundary, int lastCol, int edge, int nRows) {

    for (int count = 0, tarDiff = 3; ;
         NextDistinctPart(z, boundary, edge, tarDiff, lastCol)) {

        // index of first non-zero element
        const auto it = std::find_if(z.cbegin(), z.cend(),
                                     [](int z_i) {return z_i;});

        const int nz = std::distance(z.cbegin(), it);

        do {
            for (int k = 0; k < width; ++k) {
                mat[count + nRows * k] = z[k];
            }
            
            ++count;
        } while (std::next_permutation(z.begin() + nz, z.end()) &&
                 count < nRows);

        if (count >= nRows) {break;}
    }
}

// !mIsNull || !IncludeZero
void PartsLenPermDistinct(int* mat, std::vector<int> &z,
                          int width, int boundary, int lastCol,
                          int edge, int nRows) {
    
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

    for (int count = 0, tarDiff = 3; ;
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

template void PartsGenPermDistinct(std::vector<int>&,
                                   const std::vector<int>&,
                                   std::vector<int>&, int, int, int, int);

template void PartsGenPermDistinct(std::vector<double>&,
                                   const std::vector<double>&,
                                   std::vector<int>&, int, int, int, int);
