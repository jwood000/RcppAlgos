#include "Partitions/PartitionsMultiset.h"
#include "Partitions/PartitionsDistinct.h"
#include "Partitions/PartitionsTypes.h"
#include "Partitions/PartitionsRep.h"
#include "RMatrix.h"

template <typename T>
void PartsManager(T* mat, const std::vector<T> &v,
                  const std::vector<int> &freqs, std::vector<int> &z,
                  const PartDesign &part, int n, int m, int nRows) {
    
    int lastCol = part.width - 1;
    int boundary = lastCol;
    int edge = boundary - 1;

    switch (part.ptype) {
        case PartitionType::RepStdAll: {
            return PartsRep(mat, z, part.width, boundary, edge, lastCol, 0, nRows);
        } case PartitionType::RepNoZero: {
            return PartsRep(mat, z, part.width, boundary, edge, lastCol, 0, nRows);
        } case PartitionType::RepShort: {
            return PartsRep(mat, z, part.width, boundary, edge, lastCol, 0, nRows);
        } case PartitionType::RepCapped: {
            // return CountPartRepLenCap(part.mapTar, part.width, lenV);
        } case PartitionType::DstctStdAll: {
            // return CountPartDistinct(part.mapTar);
        } case PartitionType::DstctShort: {
            // return GetSpecialCount(part.startZ, part.mapTar, part.width);
        } case PartitionType::DstctSpecial: {
            // return GetSpecialCount(part.startZ, part.mapTar, part.width);
        } case PartitionType::DstctOneZero: {
            // return CountPartDistinctLen(part.mapTar, part.width);
        } case PartitionType::DstctNoZero: {
            // return CountPartDistinctLen(part.mapTar, part.width);
        } case PartitionType::DistCapped: {
            // return CountPartDistinctLenCap(part.mapTar, part.width, lenV);
        } case PartitionType::Multiset: {
            // return CountPartMultiset(Reps, part.startZ);
        } default: {
            // Do nothing
        }
    }
}

// template <typename T>
// void ComboParallel(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
//                    std::vector<int> &z, int n, int m, int strt, int nRows,
//                    const std::vector<int> &freqs, bool IsMult, bool IsRep) {
// 
//     if (IsMult) {
//         MultisetCombination(mat, v, z, n, m, strt, nRows, freqs);
//     } else if (IsRep) {
//         CombinationsRep(mat, v, z, n, m, strt, nRows);
//     } else {
//         CombinationsDistinct(mat, v, z, n, m, strt, nRows);
//     }
// }
// 
// template void ComboManager(int*, const std::vector<int>&,
//                            std::vector<int>&, int, int, int,
//                            const std::vector<int>&, bool, bool);
// 
// template void ComboManager(double*, const std::vector<double>&,
//                            std::vector<int>&, int, int, int,
//                            const std::vector<int>&, bool, bool);
// 
// template void ComboManager(Rbyte*, const std::vector<Rbyte>&,
//                            std::vector<int>&, int, int, int,
//                            const std::vector<int>&, bool, bool);
// 
// template void ComboManager(Rcomplex*, const std::vector<Rcomplex>&,
//                            std::vector<int>&, int, int, int,
//                            const std::vector<int>&, bool, bool);
// 
// template void ComboParallel(RcppParallel::RMatrix<int>&, const std::vector<int>&,
//                             std::vector<int>&, int, int, int, int,
//                             const std::vector<int>&, bool, bool);
// 
// template void ComboParallel(RcppParallel::RMatrix<double>&, const std::vector<double>&,
//                             std::vector<int>&, int, int, int, int,
//                             const std::vector<int>&, bool, bool);
