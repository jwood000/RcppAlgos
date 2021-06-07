#include "Partitions/PartitionsMultiset.h"
#include "Partitions/PartitionsDistinct.h"
#include "Partitions/PartitionsTypes.h"
#include "Partitions/PartitionsRep.h"
#include "RMatrix.h"

template <typename T>
void PartsGenManager(T* mat, const std::vector<T> &v,
                     std::vector<T> &partVec,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     int width, int lastCol, int lastElem, int boundary,
                     int edge, int strt, int nRows, bool IsComb,
                     bool IsRep, bool IsMult) {
    
    if (IsComb) {
        if (IsMult) {
            PartsGenMultiset(partVec, v, Reps, z, width,
                             lastElem, lastCol, nRows);
        } else if (IsRep) {
            PartsGenRep(mat, v, z, width,
                        lastElem, lastCol, strt, nRows);
        } else {
            PartsGenDistinct(mat, v, z, width,
                             lastElem, lastCol, strt, nRows);
        }
    } else {
        // if (part.ptype <= PartitionType::RepShort) {
        //     PartsRep(mat, z, part.width, boundary,
        //              edge, lastCol, strt, nRows);
        // } else if (part.ptype == PartitionType::RepCapped) {
        //     PartsGenRep(mat, v, z, part.width,
        //                 lastElem, lastCol, strt, nRows);
        // } else if (part.ptype <= PartitionType::DstctNoZero) {
        //     PartsDistinct(mat, z, part.width, boundary,
        //                   lastCol, edge, strt, nRows);
        // } else if (part.ptype == PartitionType::DistCapped) {
        //     PartsGenDistinct(mat, v, z, part.width,
        //                      lastElem, lastCol, strt, nRows);
        // } else {
        //     PartsGenMultiset(partVec, v, Reps, z, part.width,
        //                      lastElem, lastCol, nRows);
        // }
    }
}

void PartsStdManager(int* mat, std::vector<int> &z, int width,
                     int lastCol, int boundary, int edge,
                     int strt, int nRows, bool IsComb, bool IsRep) {
    
    if (IsRep && IsComb) {
        PartsRep(mat, z, width, boundary, edge, lastCol, strt, nRows);
    } else if (IsRep) {
        PartsPermRep(mat, z, width, boundary, edge, lastCol, nRows);
    } else if (IsComb) {
        PartsDistinct(mat, z, width, boundary, lastCol, edge, strt, nRows);
    } else {
        PartsPermDistinct(mat, z, width, boundary, lastCol, edge, nRows);
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

// (T* mat, const std::vector<T> &v, std::vector<T> &partVec,
//  const std::vector<int> &Reps, std::vector<int> &z,
//  const PartDesign &part, int lastCol, int lastElem,
//  int boundary, int edge, int strt, int nRows, bool IsComb)
template void PartsGenManager(int*, const std::vector<int>&,
                              std::vector<int>&, const std::vector<int>&,
                              std::vector<int>&, int, int, int, int, int,
                              int, int, bool, bool, bool);

template void PartsGenManager(double*, const std::vector<double>&,
                              std::vector<double>&, const std::vector<int>&,
                              std::vector<int>&, int, int, int, int, int,
                              int, int, bool, bool, bool);

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
