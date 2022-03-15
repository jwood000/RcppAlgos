#include "Partitions/NextPartition.h"
#include "PopulateVec.h"
#include "RMatrix.h"
#include <algorithm>  // std::next_permutation

template <typename T>
void PartsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                 int width, int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (int count = 0; count < nRows; ++count,
         NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem)) {

        for (int k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
    }
}

template <typename T>
void PartsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                 std::vector<int> &z, int strt, int width, int lastElem,
                 int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (int count = strt; count < nRows; ++count,
         NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem)) {

        for (int k = 0; k < width; ++k) {
            mat(count, k) = v[z[k]];
        }
    }
}

template <typename T>
void PartsGenPermRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                     int width, int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (int count = 0; ;
         NextRepPart(z, boundary, edge, lastCol)) {

        do {
            for (int k = 0; k < width; ++k) {
                mat[count + nRows * k] = v[z[k]];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
    }
}

template <typename T>
void PartsGenRep(std::vector<T> &partsVec, const std::vector<T> &v,
                 std::vector<int> &z, int width, int nRows, bool IsComb) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    const int lastCol = width - 1;
    const int lastElem = v.size() - 1;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (int count = 0; (edge >= 0) && (z[boundary] - z[edge] >= 2);
         NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem)) {

        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
        if (count >= nRows) break;
    }

    int count = partsVec.size() / width;

    if (count < nRows) {
        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
    }
}

void PartsRep(int* mat, std::vector<int> &z, int width,
              int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (int count = 0; count < nRows; ++count,
         NextRepPart(z, boundary, edge, lastCol)) {

        for (int k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }
    }
}

void PartsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
              int strt, int width, int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (int count = strt; count < nRows; ++count,
         NextRepPart(z, boundary, edge, lastCol)) {

        for (int k = 0; k < width; ++k) {
            mat(count, k) = z[k];
        }
    }
}

void PartsPermRep(int* mat, std::vector<int> &z, int width,
                  int lastElem, int lastCol, int nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (int count = 0; ;
         NextRepPart(z, boundary, edge, lastCol)) {

        do {
            for (int k = 0; k < width; ++k) {
                mat[count + nRows * k] = z[k];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
    }
}

template void PartsGenRep(int*, const std::vector<int>&,
                          std::vector<int>&, int, int, int, int);
template void PartsGenRep(double*, const std::vector<double>&,
                          std::vector<int>&, int, int, int, int);

template void PartsGenRep(RcppParallel::RMatrix<int>&,
                          const std::vector<int>&, std::vector<int> &z,
                          int, int, int, int, int);
template void PartsGenRep(RcppParallel::RMatrix<double>&,
                          const std::vector<double>&, std::vector<int> &z,
                          int, int, int, int, int);

template void PartsGenRep(std::vector<int>&, const std::vector<int>&,
                          std::vector<int>&, int, int, bool);
template void PartsGenRep(std::vector<double>&, const std::vector<double>&,
                          std::vector<int>&, int, int, bool);

template void PartsGenPermRep(int*, const std::vector<int>&,
                              std::vector<int>&, int, int, int, int);
template void PartsGenPermRep(double*, const std::vector<double>&,
                              std::vector<int>&, int, int, int, int);
