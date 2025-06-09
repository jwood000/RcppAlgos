#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/CompositionsDistinct.h"
#include "Partitions/PartitionsMultiset.h"
#include "Partitions/PartitionsDistinct.h"
#include "Partitions/PartitionsTypes.h"
#include "Partitions/CompositionsRep.h"
#include "Partitions/PartitionsRep.h"
#include <algorithm> // std::find
#include <numeric>
#include "RMatrix.h"

void PartsStdManager(int* mat, std::vector<int> &z, int width,
                     int lastElem, int lastCol, int nRows, bool IsComb,
                     bool IsRep, bool IsComp, bool zero_spesh) {

    // switch (ptype) {. PartitionType ptype
    //     case PartitionType::LengthOne: {
    //         if (nRows) mat[0] = z.front();
    //         break;
    //     } case PartitionType::RepStdAll: {
    //         PartsRep(mat, z, width, lastElem, lastCol, nRows);
    //         break;
    //     } case PartitionType::RepNoZero: {
    //         PartsRep(mat, z, width, lastElem, lastCol, nRows);
    //         break;
    //     } case PartitionType::RepShort: {
    //         PartsRep(mat, z, width, lastElem, lastCol, nRows);
    //         break;
    //     } case PartitionType::RepCapped: {
    //         PartsRep(mat, z, width, lastElem, lastCol, nRows);
    //         break;
    //     } case PartitionType::DstctStdAll: {
    //
    //         break;
    //     } case PartitionType::DstctMultiZero: {
    //
    //         break;
    //     } case PartitionType::DstctOneZero: {
    //
    //         break;
    //     } case PartitionType::DstctNoZero: {
    //
    //         break;
    //     } case PartitionType::DstctCapped: {
    //
    //         break;
    //     } case PartitionType::DstctCappedMZ: {
    //
    //         break;
    //     } case PartitionType::Multiset: {
    //
    //         break;
    //     } case PartitionType::CoarseGrained: {
    //
    //         break;
    //     } case PartitionType::CompRepNoZero: {
    //
    //         break;
    //     } case PartitionType::CompRepWeak: {
    //
    //         break;
    //     } case PartitionType::CmpRpZroNotWk: {
    //
    //         break;
    //     } case PartitionType::CmpDstctNoZero: {
    //         CompsRep<0>(mat, z, width, nRows);
    //         break;
    //     } case PartitionType::CmpDstctZNotWk: {
    //         CompsRep<1>(mat, z, width, nRows);
    //         break;
    //     } case PartitionType::CmpDstctMZWeak: {
    //         CompsRep<0>(mat, z, width, nRows);
    //         break;
    //     } default: {
    //
    //     }
    // }

    if (width == 1) {
        if (nRows) mat[0] = z.front();
    } else if (IsRep && IsComb) {
        PartsRep(mat, z, width, lastElem, lastCol, nRows);
    } else if (IsRep && IsComp && zero_spesh) {
        CompsRep<1>(mat, z, width, nRows);
    } else if (IsRep && IsComp) {
        CompsRep<0>(mat, z, width, nRows);
    } else if (IsComp) {
        std::vector<int> myRange(z.back());
        std::iota(myRange.begin(), myRange.end(), 1);

        std::vector<int> complement;
        std::set_difference(myRange.begin(), myRange.end(), z.begin(), z.end(),
                            std::inserter(complement, complement.begin()));

        int lastIdx = complement.size() - 1;
        const int target = std::accumulate(z.cbegin(), z.cend(), 0);
        const int numZeros = std::count(z.cbegin(), z.cend(), 0);
        int strt = 0;

        for (int i = width - numZeros, j = numZeros,
             nextStep = 0; i < width; ++i, --j) {
            nextStep += CountCompsDistinctLen(target, i);
            CompsDistinct(mat, z, complement, 0, lastIdx, z.back(),
                          target, strt, width, nextStep, nRows);
            strt = nextStep;

            for (int k = j - 1, val = 1; k < width; ++k, ++val) {
                z[k] = val;
            }

            z.back() = target - static_cast<int>((i * (i + 1)) / 2);
            complement.clear();
            std::set_difference(
                myRange.begin(), myRange.end(), z.begin(), z.end(),
                std::inserter(complement, complement.begin())
            );
        }

        CompsDistinct(mat, z, complement, 0, lastIdx, z.back(),
                      target, strt, width, nRows, nRows);
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
                     bool IsComb, bool IsRep, bool IsComp, bool zero_spesh) {

    if (width == 1) {
        if (nRows) mat[0] = v[z.front()];
    } else if (IsComb && IsRep) {
        PartsGenRep(mat, v, z, width, lastElem, lastCol, nRows);
    } else if (IsComb) {
        PartsGenDistinct(mat, v, z, width, lastElem, lastCol, nRows);
    } else if (IsRep && IsComp && zero_spesh) {
        CompsGenRep<1>(mat, v, z, width, nRows);
    } else if (IsRep && IsComp) {
        CompsGenRep<0>(mat, v, z, width, nRows);
    } else if (IsRep) {
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
                      int nRows, bool IsRep, bool IsComp, bool zero_spesh) {

    if (IsRep && IsComp && zero_spesh) {
        CompsRep<1>(mat, z, strt, width, nRows);
    } else if (IsRep && IsComp) {
        CompsRep<0>(mat, z, strt, width, nRows);
    } else if (IsRep) {
        PartsRep(mat, z, strt, width, lastElem, lastCol, nRows);
    } else {
        PartsDistinct(mat, z, strt, width, lastElem, lastCol, nRows);
    }
}

template <typename T>
void PartsGenParallel(RcppParallel::RMatrix<T> &mat,
                      const std::vector<T> &v, std::vector<int> &z,
                      int strt, int width, int lastElem, int lastCol,
                      int nRows, bool IsRep, bool IsComp, bool zero_spesh) {

    if (IsRep && IsComp && zero_spesh) {
        CompsGenRep<1>(mat, v, z, strt, width, nRows);
    } else if (IsRep && IsComp) {
        CompsGenRep<0>(mat, v, z, strt, width, nRows);
    } else if (IsRep) {
        PartsGenRep(mat, v, z, strt, width, lastElem, lastCol, nRows);
    } else {
        PartsGenDistinct(mat, v, z, strt, width, lastElem, lastCol, nRows);
    }
}

template void PartsGenManager(int*, const std::vector<int>&,
                              std::vector<int>&, int, int, int,
                              int, bool, bool, bool, bool);
template void PartsGenManager(double*, const std::vector<double>&,
                              std::vector<int>&, int, int, int,
                              int, bool, bool, bool, bool);

template void PartsGenManager(std::vector<int>&, const std::vector<int>&,
                              const std::vector<int>&, std::vector<int>&,
                              PartitionType, int, int, bool);
template void PartsGenManager(std::vector<double>&, const std::vector<double>&,
                              const std::vector<int>&, std::vector<int>&,
                              PartitionType, int, int, bool);

template void PartsGenParallel(RcppParallel::RMatrix<int>&,
                               const std::vector<int>&, std::vector<int>&,
                               int, int, int, int, int, bool, bool, bool);
template void PartsGenParallel(RcppParallel::RMatrix<double>&,
                               const std::vector<double>&, std::vector<int>&,
                               int, int, int, int, int, bool, bool, bool);
