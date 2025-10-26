#include "Partitions/NextPartition.h"
#include "PopulateVec.h"
#include "RMatrix.h"
#include <algorithm>  // std::next_permutation

template <typename T>
int PartsGenRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                std::size_t width, int lastElem,
                int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (std::size_t count = 0; count < nRows; ++count,
         NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = v[z[k]];
        }
    }

    return 1;
}

template <typename T>
int PartsGenRep(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
                std::vector<int> &z, int strt, std::size_t width,
                int lastElem, int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (std::size_t count = strt; count < nRows; ++count,
         NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = v[z[k]];
        }
    }

    return 1;
}

template <typename T>
int PartsGenPermRep(T* mat, const std::vector<T> &v, std::vector<int> &z,
                    std::size_t width, int lastElem,
                    int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (std::size_t count = 0; ;
         NextRepPart(z, boundary, edge, lastCol)) {

        do {
            for (std::size_t k = 0; k < width; ++k) {
                mat[count + nRows * k] = v[z[k]];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
    }

    return 1;
}

template <typename T>
int PartsGenRep(std::vector<T> &partsVec, const std::vector<T> &v,
                std::vector<int> &z, std::size_t width,
                std::size_t nRows, bool IsComb) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    const int lastCol = width - 1;
    const int lastElem = v.size() - 1;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (std::size_t count = 0; (edge >= 0) && (z[boundary] - z[edge] >= 2);
         NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem)) {

        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
        if (count >= nRows) break;
    }

    std::size_t count = partsVec.size() / width;

    if (count < nRows) {
        PopulateVec(v, partsVec, z, count, width, nRows, IsComb);
    }

    return 1;
}

int PartsRep(int* mat, std::vector<int> &z, std::size_t width,
             int lastElem, int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (std::size_t count = 0; count < nRows; ++count,
         NextRepPart(z, boundary, edge, lastCol)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat[count + nRows * k] = z[k];
        }
    }

    return 1;
}

int PartsRep(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
             int strt, std::size_t width, int lastElem,
             int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (std::size_t count = strt; count < nRows; ++count,
         NextRepPart(z, boundary, edge, lastCol)) {

        for (std::size_t k = 0; k < width; ++k) {
            mat(count, k) = z[k];
        }
    }

    return 1;
}

int PartsPermRep(int* mat, std::vector<int> &z, std::size_t width,
                 int lastElem, int lastCol, std::size_t nRows) {

    int edge = 0;
    int pivot = 0;
    int boundary = 0;

    PrepareRepPart(z, boundary, pivot, edge, lastElem, lastCol);

    for (std::size_t count = 0; ;
         NextRepPart(z, boundary, edge, lastCol)) {

        do {
            for (std::size_t k = 0; k < width; ++k) {
                mat[count + nRows * k] = z[k];
            }

            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < nRows);

        if (count >= nRows) {break;}
    }

    return 1;
}

template int PartsGenRep(int*, const std::vector<int>&,
                         std::vector<int>&, std::size_t,
                         int, int, std::size_t);
template int PartsGenRep(double*, const std::vector<double>&,
                         std::vector<int>&, std::size_t,
                         int, int, std::size_t);

template int PartsGenRep(RcppParallel::RMatrix<int>&,
                         const std::vector<int>&, std::vector<int> &z,
                         int, std::size_t, int, int, std::size_t);
template int PartsGenRep(RcppParallel::RMatrix<double>&,
                         const std::vector<double>&, std::vector<int> &z,
                         int, std::size_t, int, int, std::size_t);

template int PartsGenRep(std::vector<int>&, const std::vector<int>&,
                         std::vector<int>&, std::size_t, std::size_t, bool);
template int PartsGenRep(std::vector<double>&, const std::vector<double>&,
                         std::vector<int>&, std::size_t, std::size_t, bool);

template int PartsGenPermRep(int*, const std::vector<int>&,
                             std::vector<int>&, std::size_t,
                             int, int, std::size_t);
template int PartsGenPermRep(double*, const std::vector<double>&,
                             std::vector<int>&, std::size_t,
                             int, int, std::size_t);
