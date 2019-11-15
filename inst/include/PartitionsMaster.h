#ifndef PARTITIONS_MASTER_H
#define PARTITIONS_MASTER_H

#include "GeneralPartitions.h"
#include "Cpp14MakeUnique.h"

namespace Partitions {
    
    void NextDistinct(std::vector<int> &z, int &boundary,
                      int &edge, int &tarDiff, int lastCol) {
        
        if (z[boundary] - z[edge] != tarDiff)
            boundary = edge + 1;
        
        ++z[edge];
        --z[boundary];
        int myEdgePlus = z[edge] + boundary - edge;
        
        for (; boundary < lastCol; ++boundary, ++myEdgePlus) {
            z[lastCol] += (z[boundary] - myEdgePlus);
            z[boundary] = myEdgePlus;
        }
        
        while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2) {
            --boundary;
        }
        
        edge = boundary - 1;
        tarDiff = 3;
        
        while (edge && (z[boundary] - z[edge]) < tarDiff) {
            --edge;
            ++tarDiff;
        }
    }
    
    template <typename typeRcpp>
    void PartsStdDistinct(typeRcpp &partitionsMatrix, std::vector<int> &z, std::size_t m,
                          std::size_t numRows, bool getAll, int boundary, int lastCol,
                          int edge, bool IsComb, bool IncludeZero, bool mIsNull) {
        int tarDiff = 3;
        
        if (IsComb) {
            for (std::size_t j = 0, lim = numRows - 1; j < lim; ++j) {
                for (std::size_t k = 0; k < m; ++k)
                    partitionsMatrix(j, k) = z[k];
                
                NextDistinct(z, boundary, edge, tarDiff, lastCol);
            }
            
            for (std::size_t k = 0; k < m; ++k)
                partitionsMatrix(numRows - 1, k) = z[k];
        } else {
            if (mIsNull && IncludeZero) {
                for (std::size_t count = 0; ; NextDistinct(z, boundary, edge, tarDiff, lastCol)) {
                    
                    const auto it = std::find_if(z.cbegin(), z.cend(),
                                                 [](int z_i) {return z_i;});
                    const int firstNonZero = std::distance(z.cbegin(), it);

                    do {
                        for (std::size_t k = 0; k < m; ++k)
                            partitionsMatrix(count, k) = z[k];

                        ++count;
                    } while (std::next_permutation(z.begin() + firstNonZero, z.end()) && count < numRows);

                    if (count >= numRows) {break;}
                }
            } else {
                const std::size_t indexRows = static_cast<std::size_t>(NumPermsNoRep(m, m));
                auto indexMat = FromCpp14::make_unique<int[]>(indexRows * m);
                
                std::vector<int> indexVec(m);
                std::iota(indexVec.begin(), indexVec.end(), 0);
                
                for (std::size_t i = 0, myRow = 0; i < indexRows; ++i, myRow += m) {
                    for (std::size_t j = 0; j < m; ++j)
                        indexMat[myRow + j] = indexVec[j];
                    
                    std::next_permutation(indexVec.begin(), indexVec.end());
                }
                
                for (std::size_t count = 0; ; NextDistinct(z, boundary, edge, tarDiff, lastCol)) {
                    for (std::size_t j = 0, myRow = 0; j < indexRows; ++count, ++j)
                        for (std::size_t k = 0; k < m; ++k, ++myRow)
                            partitionsMatrix(count, k) = z[indexMat[myRow]];
                    
                    if (count >= numRows) {break;}
                }
            }
        }
    }

    void NextPartition(std::vector<int> &z, int &boundary, int &edge, int lastCol) {
        
        if (z[boundary] - z[edge] != 2)
            boundary = edge + 1;
        
        ++z[edge];
        --z[boundary];
        const int myEdge = z[edge];
        
        for (; boundary < lastCol; ++boundary) {
            z[lastCol] += (z[boundary] - myEdge);
            z[boundary] = myEdge;
        }
        
        const int currMax = z[boundary];
        
        while (boundary > 1 && z[boundary - 1] == currMax)
            --boundary;
        
        edge = boundary - 1;
        const int edgeTest = z[boundary] - 2;
        
        while (edge && edgeTest < z[edge])
            --edge;
    }
    
    template <typename typeRcpp>
    void PartitionsStandard(typeRcpp &partitionsMatrix, std::vector<int> &z, std::size_t m,
                            std::size_t numRows, bool IsComb, int boundary, int edge, 
                            int lastCol, bool includeZero, bool mIsNull) {
        if (IsComb) {
            if (includeZero) {
                for (std::size_t j = 0, lim = numRows - 1; j < lim; ++j, NextPartition(z, boundary, edge, lastCol))
                    for (int k = lastCol; k >= 0 && z[k]; --k)
                        partitionsMatrix(j, k) = z[k];
            } else {
                for (std::size_t j = 0, lim = numRows - 1; j < lim; ++j, NextPartition(z, boundary, edge, lastCol))
                    for (std::size_t k = 0; k < m; ++k)
                        partitionsMatrix(j, k) = z[k];
            }
            
            for (std::size_t k = 0; k < m; ++k)
                partitionsMatrix(numRows - 1, k) = z[k];
        } else {
            if (mIsNull && includeZero) {
                for (std::size_t count = 0;; NextPartition(z, boundary, edge, lastCol)) {
                    
                    const auto it = std::find_if(z.cbegin(), z.cend(), 
                                                 [](int z_i) {return z_i;});
                    const int firstNonZero = std::distance(z.cbegin(), it);
                    
                    do {
                        for (int k = lastCol; k >= 0 && z[k]; --k)
                            partitionsMatrix(count, k) = z[k];
                        
                        ++count;
                    } while (std::next_permutation(z.begin() + firstNonZero, z.end()) && count < numRows);
                    
                    if (count >= numRows) {break;}
                }
            } else {
                for (std::size_t count = 0;; NextPartition(z, boundary, edge, lastCol)) {
                    do {
                        for (std::size_t k = 0; k < m; ++k)
                            partitionsMatrix(count, k) = z[k];
                        
                        ++count;
                    } while (std::next_permutation(z.begin(), z.end()) && count < numRows);
                    
                    if (count >= numRows) {break;}
                }
            }
        }
    }
        
    template <typename typeRcpp>
    typeRcpp GeneralPartitions(const std::vector<int64_t> &v, std::vector<int> &z, const std::vector<int> &Reps,
                               PartitionType PartType, int64_t target, int lenV, int m, bool IsRep, bool IsMult,
                               double numParts, bool IsComb, bool xtraCol, bool bUserRows, bool mIsNull, bool getAll) {

        const bool IncludeZero = (v.front() == 0);
        
        // When Rm is NULL, m = n, and we have zero, then we need
        // to reduce m by 1. In this case lenV = target + 1
        if (mIsNull && IncludeZero && m == lenV)
            m = lenV - 1;
        
        std::size_t nCols = (xtraCol) ? m + 1 : m;
        
        const int lastElem = lenV - 1;
        int lastCol = m - 1;
        
        if (PartType > PartitionType::PartGeneral) {
            const bool IsDistinct = (PartType > PartitionType::PartTradNoZero);
            int boundary = lastCol;
            int edge = boundary - 1;
            
            // If we have determined that the Partitions are of standard
            // form, we have already checked that the number of results
            // is within integer range in the calling code.
            const int nRows = static_cast<int>(numParts);
            nCols = (xtraCol) ? m + 1 : m;
            
            typeRcpp partitionsMatrix = (IncludeZero && (IsComb || mIsNull)) ? 
                                            typeRcpp(nRows, nCols) : 
                                                Rcpp::no_init_matrix(nRows, nCols);
            
            if (IsDistinct) {
                PartsStdDistinct(partitionsMatrix, z, m, nRows, getAll,
                                 boundary, lastCol, edge, IsComb, IncludeZero, mIsNull);
            } else {
                PartitionsStandard(partitionsMatrix, z, m, nRows, IsComb,
                                   boundary, edge, lastCol, IncludeZero, mIsNull);
            }
            
            if (xtraCol)
                for (int i = 0; i < nRows; ++i)
                    partitionsMatrix(i, m) = target;
            
            return partitionsMatrix;
        }
        
        std::vector<int64_t> partitionsVec;
        const double vecMax = partitionsVec.max_size() / m;
        const double upperBound = std::min(vecMax, static_cast<double>(std::numeric_limits<int>::max()));
        const int maxRows = std::min(upperBound, static_cast<double>(numParts));
        
        if (bUserRows)
            partitionsVec.reserve(maxRows * m);
        
        int result = 0;
        
        if (IsRep) {
            result = PartitionsRep(m, v, target, lastElem, lastCol, maxRows, IsComb, partitionsVec);
        } else if (IsMult) {
            result = PartitionsMultiSet(m, v, target, lastElem, lastCol, maxRows, IsComb, Reps, partitionsVec);
        } else {
            result = PartitionsDistinct(m, v, target, lastElem, lastCol, maxRows, IsComb, partitionsVec);
        }
        
        if (result) {
            std::size_t partitionLen = partitionsVec.size();
            std::size_t numResult = partitionLen / m;
            typeRcpp partitionsMatrix = Rcpp::no_init_matrix(numResult, nCols);
            
            for (std::size_t i = 0, k = 0; i < numResult; ++i)
                for (int j = 0; j < m; ++j, ++k)
                    partitionsMatrix(i, j) = partitionsVec[k];
                
            if (xtraCol)
                for (std::size_t i = 0; i < numResult; ++i)
                    partitionsMatrix(i, m) = target;
                    
            if (partitionLen >= upperBound) {
                Rcpp::warning("The algorithm terminated early as the number of results "
                              "meeting the criteria exceeds the container's maximum "
                              "capacity or 2^31 - 1");
            }
            
            return partitionsMatrix;
        } else {
            typeRcpp trivialRet(0, m);
            return trivialRet;
        }
    }
}

#endif
