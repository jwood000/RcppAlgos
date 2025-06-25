#include "cpp11/protect.hpp"

#include "Partitions/CompositionsDistinct.h"
#include "Partitions/PartitionsMultiset.h"
#include "Partitions/PartitionsDistinct.h"
#include "Partitions/PartitionsTypes.h"
#include "Partitions/CompositionsRep.h"
#include "Partitions/PartitionsRep.h"
#include <algorithm> // std::find
#include <numeric>
#include "RMatrix.h"

int PartsStdManager(int* mat, std::vector<int> &z, int width, int lastElem,
                    int lastCol, int nRows, PartitionType ptype) {

    switch (ptype) {
        case PartitionType::LengthOne: {
            if (nRows) mat[0] = z.front();
            return 1;
        } case PartitionType::RepStdAll: {
            PartsRep(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::RepNoZero: {
            PartsRep(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::RepShort: {
            PartsRep(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::DstctStdAll: {
            PartsDistinct(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::DstctMultiZero: {
            PartsDistinct(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::DstctOneZero: {
            PartsDistinct(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::DstctNoZero: {
            PartsDistinct(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::CompRepNoZero: {
            CompsRep<1>(mat, z, width, nRows);
            return 1;
        } case PartitionType::CompRepWeak: {
            CompsRep<0>(mat, z, width, nRows);
            return 1;
        } case PartitionType::CmpRpZroNotWk: {
            CompsRep<1>(mat, z, width, nRows);
            return 1;
        } case PartitionType::PrmRepPartNoZ: {
            PartsPermRep(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::PrmRepPart: {
            PartsPermRep(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::PrmDstPartNoZ: {
            PartsPermDistinct(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::PrmDstPrtOneZ: {
            PartsPermDistinct(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::PrmDstPartMZ: {
            PartsPermZeroDistinct(mat, z, width, lastElem, lastCol, nRows);
            return 1;
        } case PartitionType::CmpDstctNoZero: {
            CompsDistinct(mat, z, width, nRows);
            return 1;
        } case PartitionType::CmpDstctZNotWk: {
            CompsDistinct(mat, z, width, nRows);
            return 1;
        } case PartitionType::CmpDstctMZWeak: {
            CompsDistinct(mat, z, width, nRows);
            return 1;
        } default: {
            // This should not happen
            cpp11::stop("Casse not supported");
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
                     PartitionType ptype, int width, int nRows) {

    if (width == 1) {
        if (nRows) partsVec.push_back(v[z.front()]);
    } else if (ptype == PartitionType::Multiset) {
        PartsGenMultiset(partsVec, v, Reps, z, width, nRows, true);
    } else if (ptype == PartitionType::PrmMultiset) {
        PartsGenMultiset(partsVec, v, Reps, z, width, nRows, false);
    } else if (ptype == PartitionType::RepCapped) {
        PartsGenRep(partsVec, v, z, width, nRows, true);
    } else if (ptype == PartitionType::PrmRepCapped) {
        PartsGenRep(partsVec, v, z, width, nRows, false);
    } else if (ptype == PartitionType::DstctCapped) {
        PartsGenDistinct(partsVec, v, z, width, nRows, true);
    } else {
        PartsGenDistinct(partsVec, v, z, width, nRows, false);
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
                              PartitionType, int, int);
template void PartsGenManager(std::vector<double>&, const std::vector<double>&,
                              const std::vector<int>&, std::vector<int>&,
                              PartitionType, int, int);

template void PartsGenParallel(RcppParallel::RMatrix<int>&,
                               const std::vector<int>&, std::vector<int>&,
                               int, int, int, int, int, bool, bool, bool);
template void PartsGenParallel(RcppParallel::RMatrix<double>&,
                               const std::vector<double>&, std::vector<int>&,
                               int, int, int, int, int, bool, bool, bool);
