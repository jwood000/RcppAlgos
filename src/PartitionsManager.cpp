#include "Partitions/PartitionsMultiset.h"
#include "Partitions/PartitionsDistinct.h"
#include "Partitions/PartitionsTypes.h"
#include "Partitions/PartitionsRep.h"
#include <algorithm> // std::find
#include "RMatrix.h"

void PartsStdManager(int* mat, std::vector<int> &z, int width,
                     int lastElem, int lastCol, int nRows,
                     bool IsComb, bool IsRep) {

    if (width == 1) {
        if (nRows) mat[0] = z.front();
    } else if (IsRep && IsComb) {
        PartsRep(mat, z, width, lastElem, lastCol, nRows);
    } else if (IsRep) {
        PartsPermRep(mat, z, width, lastElem, lastCol, nRows);
    } else if (IsComb) {
        PartsDistinct(mat, z, width, lastElem, lastCol, nRows);
    } else {
        const auto it = std::find(z.rbegin(), z.rend(), 0);
        const int dist = std::distance(it, z.rend());

        if (dist > 1) {
            PartsPermZeroDistinct(mat, z, width, lastElem, lastCol, nRows);
        } else {
            PartsPermDistinct(mat, z, width, lastElem, lastCol, nRows);
        }
    }
}

template <typename T>
void PartsGenManager(T* mat, const std::vector<T> &v, std::vector<int> &z,
                     int width, int lastElem, int lastCol, int nRows,
                     bool IsComb, bool IsRep) {

    if (width == 1) {
        if (nRows) mat[0] = v[z.front()];
    } else if (IsComb) {
        if (IsRep) {
            PartsGenRep(mat, v, z, width, lastElem, lastCol, nRows);
        } else {
            PartsGenDistinct(mat, v, z, width, lastElem, lastCol, nRows);
        }
    } else {
        if (IsRep) {
            PartsGenPermRep(mat, v, z, width, lastElem, lastCol, nRows);
        } else {
            const auto it = std::find(z.rbegin(), z.rend(), 0);
            const int dist = std::distance(it, z.rend());

            if (dist > 1) {
                PartsGenPermZeroDistinct(mat, v, z, width,
                                         lastElem, lastCol, nRows);
            } else {
                PartsGenPermDistinct(mat, v, z, width,
                                     lastElem, lastCol, nRows);
            }
        }
    }
}

template <typename T>
void PartsGenManager(std::vector<T> &partsVec, const std::vector<T> &v,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     PartitionType ptype, int width, int nRows,
                     bool IsComb) {

    if (width == 1) {
        if (nRows) partsVec.push_back(v[z.front()]);
    } else if (ptype == PartitionType::Multiset) {
        PartsGenMultiset(partsVec, v, Reps, z, width, nRows, IsComb);
    } else if (ptype == PartitionType::RepCapped) {
        PartsGenRep(partsVec, v, z, width, nRows, IsComb);
    } else {
        PartsGenDistinct(partsVec, v, z, width, nRows, IsComb);
    }
}

void PartsStdParallel(RcppParallel::RMatrix<int> &mat, std::vector<int> &z,
                      int strt, int width, int lastElem, int lastCol,
                      int nRows, bool IsRep) {

    if (IsRep) {
        PartsRep(mat, z, strt, width, lastElem, lastCol, nRows);
    } else {
        PartsDistinct(mat, z, strt, width, lastElem, lastCol, nRows);
    }
}

template <typename T>
void PartsGenParallel(RcppParallel::RMatrix<T> &mat,
                      const std::vector<T> &v, std::vector<int> &z,
                      int strt, int width, int lastElem, int lastCol,
                      int nRows, bool IsRep) {

    if (IsRep) {
        PartsGenRep(mat, v, z, strt, width, lastElem, lastCol, nRows);
    } else {
        PartsGenDistinct(mat, v, z, strt, width, lastElem, lastCol, nRows);
    }
}

template void PartsGenManager(int*, const std::vector<int>&,
                              std::vector<int>&, int, int,
                              int, int, bool, bool);
template void PartsGenManager(double*, const std::vector<double>&,
                              std::vector<int>&, int, int,
                              int, int, bool, bool);

template void PartsGenManager(std::vector<int>&, const std::vector<int>&,
                              const std::vector<int>&, std::vector<int>&,
                              PartitionType, int, int, bool);
template void PartsGenManager(std::vector<double>&, const std::vector<double>&,
                              const std::vector<int>&, std::vector<int>&,
                              PartitionType, int, int, bool);

template void PartsGenParallel(RcppParallel::RMatrix<int>&,
                               const std::vector<int>&, std::vector<int>&,
                               int, int, int, int, int, bool);
template void PartsGenParallel(RcppParallel::RMatrix<double>&,
                               const std::vector<double>&, std::vector<int>&,
                               int, int, int, int, int, bool);
