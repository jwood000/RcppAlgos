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

    void DistinctIntermediates(std::size_t r, double numRows, std::vector<int> &z,
                               int boundary, int edge, int lastCol, bool CombOnly,
                               bool mIsNull, std::vector<int64_t> &partitionsVec) {
        int tarDiff = 3;
        std::size_t count = 0;
        
        if (CombOnly) {
            for (; edge >= 0 && (z[boundary] - z[edge]) >= tarDiff && count < numRows; ++count) {
                for (std::size_t k = 0; k < r; ++k)
                    partitionsVec.push_back(z[k]);
                
                NextDistinct(z, boundary, edge, tarDiff, lastCol);
            }
        } else {
            if (mIsNull) {
                for (; edge >= 0 && (z[boundary] - z[edge]) >= tarDiff; 
                            NextDistinct(z, boundary, edge, tarDiff, lastCol)) {
                    
                    auto it = std::find_if(z.cbegin(), z.cend(),
                                           [](int z_i) {return z_i;});
                    const int firstNonZero = std::distance(z.cbegin(), it);
                    
                    do {
                        for (std::size_t k = 0; k < r; ++k)
                            partitionsVec.push_back(z[k]);
                        
                        ++count;
                    } while (std::next_permutation(z.begin() + firstNonZero, z.end()) && count < numRows);
                    
                    if (count >= numRows) {break;}
                }
            } else {
                for (; edge >= 0 && (z[boundary] - z[edge]) >= tarDiff; 
                            NextDistinct(z, boundary, edge, tarDiff, lastCol)) {
                    do {
                        for (std::size_t k = 0; k < r; ++k)
                            partitionsVec.push_back(z[k]);
                        
                        ++count;
                    } while (std::next_permutation(z.begin(), z.end()) && count < numRows);
                    
                    if (count >= numRows) {break;}
                }
            }
        }
        
        if (count < numRows) {
            if (CombOnly) {
                for (std::size_t k = 0; k < r; ++k)
                    partitionsVec.push_back(z[k]);
            } else {
                do {
                    for (std::size_t k = 0; k < r; ++k)
                        partitionsVec.push_back(z[k]);
                    
                    ++count;
                } while (std::next_permutation(z.begin(), z.end()) && count < numRows);
            }
        }
    }
    
    template <typename typeRcpp>
    void PartsStdDistinct(std::size_t r, std::size_t numRows, typeRcpp &partitionsMatrix, bool getAll,
                          int boundary, int edge, int lastCol, std::vector<int> &z,
                          bool isComb, bool includeZero, bool mIsNull, const std::vector<int64_t> &partitionsVec) {
        if (isComb) {
            if (getAll) {
                int tarDiff = 3;
                
                for (std::size_t j = 0, lim = numRows - 1; j < lim; ++j) {
                    for (std::size_t k = 0; k < r; ++k)
                        partitionsMatrix(j, k) = z[k];
                    
                    NextDistinct(z, boundary, edge, tarDiff, lastCol);
                }
                
                for (std::size_t k = 0; k < r; ++k)
                    partitionsMatrix(numRows - 1, k) = z[k];
            } else {
                for (std::size_t i = 0, k = 0; i < numRows; ++i)
                    for (std::size_t j = 0; j < r; ++j, ++k)
                        partitionsMatrix(i, j) = partitionsVec[k];
            }
        } else {
            if (includeZero) {
                for (std::size_t i = 0, k = 0; i < numRows; ++i)
                    for (std::size_t j = 0; j < r; ++j, ++k)
                        partitionsMatrix(i, j) = partitionsVec[k];
            } else {
                const std::size_t partSize = partitionsVec.size() / r;
                std::size_t indexRows = static_cast<std::size_t>(NumPermsNoRep(r, r));
                auto indexMat = FromCpp14::make_unique<int[]>(indexRows * r);
                
                std::vector<int> indexVec(r);
                std::iota(indexVec.begin(), indexVec.end(), 0);
                
                for (std::size_t i = 0, myRow = 0; i < indexRows; ++i, myRow += r) {
                    for (std::size_t j = 0; j < r; ++j)
                        indexMat[myRow + j] = indexVec[j];
                    
                    std::next_permutation(indexVec.begin(), indexVec.end());
                }
                
                for (std::size_t i = 0, count = 0; i < partSize; ++i)
                    for (std::size_t j = 0, myRow = 0, base = i * r; j < indexRows; ++count, ++j)
                        for (std::size_t k = 0; k < r; ++k, ++myRow)
                            partitionsMatrix(count, k) = partitionsVec[base + indexMat[myRow]];
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
    void PartitionsStandard(std::size_t r, std::size_t numRows, typeRcpp &partitionsMatrix,
                            std::vector<int> &z, bool isComb, int boundary, int edge, 
                            int lastCol, bool includeZero, bool mIsNull) {
        if (isComb) {
            if (includeZero) {
                for (std::size_t j = 0, lim = numRows - 1; j < lim; ++j, NextPartition(z, boundary, edge, lastCol))
                    for (int k = lastCol; k >= 0 && z[k]; --k)
                        partitionsMatrix(j, k) = z[k];
            } else {
                for (std::size_t j = 0, lim = numRows - 1; j < lim; ++j, NextPartition(z, boundary, edge, lastCol))
                    for (std::size_t k = 0; k < r; ++k)
                        partitionsMatrix(j, k) = z[k];
            }
            
            for (std::size_t k = 0; k < r; ++k)
                partitionsMatrix(numRows - 1, k) = z[k];
        } else {
            if (mIsNull && includeZero) {
                for (std::size_t count = 0;; NextPartition(z, boundary, edge, lastCol)) {
                    
                    auto it = std::find_if(z.cbegin(), z.cend(), 
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
                        for (std::size_t k = 0; k < r; ++k)
                            partitionsMatrix(count, k) = z[k];
                        
                        ++count;
                    } while (std::next_permutation(z.begin(), z.end()) && count < numRows);
                    
                    if (count >= numRows) {break;}
                }
            }
        }
    }
    
    // Credit to Robin K. S. Hankin, author of the excellent partitions package.
    // From the partitions.c, here are Hankin's comments for c_numbdiffparts:
    //      "the recursion on p826 of Abramowitz and Stegun"
    double c_numbdiffparts(int n) {
        std::vector<double> qq(n, 1);
        
        for(int i = 2 ; i < n; ++i) {
            qq[i] = 0;
            
            for (int s = 1, f = 5, r = 2; i >= r; r += f, f += 3, s *= -1) {
                qq[i] += s * qq[i - r];
                if(i == r * 2) {qq[i] -= s;}
            }
            
            for (int s = 1, f = 4, r = 1; i >= r; r += f, f += 3, s *= -1) {
                qq[i] += s * qq[i - r];
                if(i == r * 2) {qq[i] -= s;}
            }
        }
        
        return qq.back();
    }
    
    // Credit to user @m69 on stackoverflow:
    //       https://stackoverflow.com/a/32918426/4408538
    // *********************************************************************** //
    Rcpp::NumericMatrix memoize;
    
    double p(int n, int r) {
        if (n <= r + 1) return 1;
        if (memoize(n - r, r - 2)) return memoize(n - r, r - 2);
        int myLim = n / r;
        if (r == 2) return myLim;
        double count = 0;
        for (; myLim--; n -= r) count += (memoize(n - r, r - 3) = p(n - 1, r - 1));
        return count;
    }
    
    double partitionCount(int n, int r, bool includeZero) {
        double count = 0;
        
        if (includeZero) {
            // Remember, if zero is included, the vector has length
            // n + 1, therefore we must subtract 1 to count properly
            const int n1 = n - 1;
            Rcpp::NumericMatrix bigRefill(n1, r);
            memoize = bigRefill;
            
            for (int i = r; i > 1; --i) {
                for (int j = 0; j < n1 - i + 1; ++j) {
                    for (int k = 0; k < i; ++k) {
                        memoize(j, k) = 0;
                    }
                }
                
                count += p(n1, i);
            }
            
            ++count;  // Add 1 for the case p(n, 1)
        } else {
            Rcpp::NumericMatrix refill(n - r + 1, r); // Initialize matrix to zero
            memoize = refill;
            count = p(n, r);
        }
        
        return count;
    }
    
    // *********************************************************************** //
    
    struct distinctType {
        int limit = 0;
        bool getAll = false;
    };
    
    distinctType DistinctAttr(int n, int r, bool isRep, bool isMult, int64_t target,
                              const std::vector<int> &Reps, bool includeZero) {
        int limit = 0;
        bool getAll = false;
        
        if (isMult || !isRep) {
            // The eqn below can be derived by taking note that the
            // smallest number of elements whose sum is at least 
            // the target will be comprised of the first x numbers.
            // That is, we need to solve for x such that:
            //
            //        sum(1:(x - 1)) <= target <= sum(1:x)
            //
            // These are triangle numbers which have the form:
            //
            //              sum(1:x) = x * (x + 1) / 2
            //
            // Given n = target, we have:
            //
            //    x * (x + 1) / 2 >= n  -->>  x^2 + x - 2n >= 0
            //
            // Finally, using the quadratic formula, we obtain:
            //
            //      x = (-1 + sqrt(1 + 4 * 1 * 2n)) / 2 * 1
            //
            // After solving for x, if sum(1:x) > target, we know
            // that the solution with the fewest number of elements
            // will contain x - 1 elements, hence std::floor.
            
            limit = static_cast<int>(std::floor((-1 + std::sqrt(1 + 8 * n)) / 2));
            
            if (isMult) {
                // Ensure all elements except the first element are 1. The first
                // element should be zero and thus have a higher frequency in
                // order to test for partitions of different length.
                
                bool allOne = std::all_of(Reps.cbegin() + 1, Reps.cend(), 
                                          [](int v_i) {return v_i == 1;});
                
                if (includeZero && n == target + 1 && allOne) {
                    if (r >= limit) {
                        if (Reps[0] >= (limit - 1)) {getAll = true;}
                    } else {
                        limit = r;
                    }
                } else {
                    limit = 0;
                }
            } else if (!isRep) {
                limit = std::min(r, limit);
            }
        }
        
        // N.B. if limit = 0, this means we either have isRep = true,
        // or we are not going to use the optimized algorithm. In this
        // case, we revert to the general algorithm.
        
        distinctType res;
        res.limit = limit;
        res.getAll = getAll;
        
        return res;
    }
        
    template <typename typeRcpp>
    typeRcpp GeneralPartitions(int n, int r, std::vector<int64_t> &v, int64_t target, bool isRep,
                               bool isMult, std::vector<int> &Reps, double numRows, bool isComb,
                               bool xtraCol, bool bUserRows, bool mIsNull) {
        
        std::vector<int64_t> partitionsVec;
        const double vecMax = partitionsVec.max_size() / r;
        if (!bUserRows) numRows = vecMax;
        
        if (isMult) {
            for (int i = 0; i < (n - 1); ++i) {
                for (int j = (i + 1); j < n; ++j) {
                    if (v[i] > v[j]) {
                        std::swap(v[i], v[j]);
                        std::swap(Reps[i], Reps[j]);
                    }
                }
            }
        } else {
            std::sort(v.begin(), v.end());
        }
        
        const bool includeZero = (v.front() == 0);
        if (mIsNull && includeZero) {r = n - 1;}
        int64_t myMax = v.back();
        std::size_t nCols = (xtraCol) ? r + 1 : r;
        
        const int lastElem = n - 1;
        int lastCol = r - 1;
        
        // We don't need to check the edge case when n == target and there is a
        // zero. Remember, n is the length of the vector v, so we could have a
        // situation where v = c(0, 2, 3, 4, 5) -->> length(v) = 5. This would
        // cause a problem if we were to allow this through, however, in the
        // calling code (i.e. Combinatorics.cpp), we ensure that the distance
        // between each element is the same. This means for the example we gave,
        // we would have length(unique(diff(v))) > 1, which means PartitionCase
        // (see Combinatorics.cpp) would be false, and thus the general 
        // algorithm would be executed.
        //
        // We do have to ensure that the smallest element is non-negative, othe-
        // rwise, cases like v = seq(-8, 10, 2), m = 7, rep = TRUE, & limit = 10
        // would pass as v = 0:9, m = 7, rep = TRUE, & limit = 9, --or--
        // v = 1:10, m = 7, rep = TRUE, & limit = 10
        
        if (myMax == target && (n == target || n == (target + 1)) && v.front() >= 0) {
            
            distinctType distinctTest = DistinctAttr(n, r, isRep, isMult, target, Reps, includeZero);
            const bool isDistinct = distinctTest.limit > 0;
            
            if (isDistinct || isRep) {
                if (isDistinct) {lastCol = distinctTest.limit - 1;}
                const std::size_t uR = (isDistinct) ? distinctTest.limit : r;
                int boundary = lastCol;
                int edge = boundary - 1;
                std::vector<int> z(lastCol + 1, 0);
                
                if (isDistinct) {
                    if (includeZero) {
                        if (isMult) {
                            if (distinctTest.getAll) {
                                z[lastCol] = lastElem;
                            } else {
                                if (Reps[0] >= static_cast<int>(uR - 1)) {
                                    z[lastCol] = lastElem;
                                } else {
                                    std::iota(z.begin() + Reps[0] - 1, z.end(), 0);
                                    z[lastCol] = target - (uR - Reps[0]) * (uR - Reps[0] - 1) / 2;
                                }
                            }
                        } else {
                            std::iota(z.begin(), z.end(), 0);
                            z[lastCol] = target - (uR - 1) * (uR - 2) / 2;
                        }
                    } else {
                        std::iota(z.begin(), z.end(), 1);
                        z[lastCol] = target - uR * (uR - 1) / 2;
                    }
                } else {
                    if (includeZero) {
                        z[lastCol] = lastElem;
                    } else {
                        std::fill(z.begin(), z.end(), 1);
                        z[lastCol] = n - uR + 1;
                    }
                }
                
                double partCountTest = 0;
                
                if (isDistinct) {
                    if (distinctTest.getAll) {
                        partCountTest = c_numbdiffparts(n);
                    } else if (isComb || !includeZero) {  // **DistinctLenR
                        DistinctIntermediates(uR, numRows, z, boundary, edge,
                                              lastCol, true, mIsNull, partitionsVec);
                        partCountTest = partitionsVec.size() / uR;
                    }
                    
                    if (!isComb) {
                        if (includeZero) {
                            // Given z = c(0, 0, 3, 4)
                            // When m is Null, here is the output (only permute last
                            // 2 indices):
                            //                  c(0, 0, 3, 4)
                            //                  c(0, 0, 4, 3)
                            //
                            // And when m is give, we have (permute all indices):
                            //                   [,1] [,2] [,3] [,4]
                            //              [1,]    0    0    3    4
                            //              [2,]    0    0    4    3
                            //              [3,]    0    3    0    4
                            //                 .    .    .    .    .
                            //             [10,]    4    0    0    3
                            //             [11,]    4    0    3    0
                            //             [12,]    4    3    0    0
                            
                            DistinctIntermediates(uR, numRows, z, boundary, edge,
                                                  lastCol, false, mIsNull, partitionsVec);
                            partCountTest = partitionsVec.size() / uR;
                        } else {
                            // In this case we have already determined the number
                            // of partitions of length r above (see DistinctLenR).
                            partCountTest *= NumPermsNoRep(r, r);
                        }
                    }
                } else {
                    if (isComb) {
                        partCountTest = partitionCount(n, r, includeZero);
                    } else if (mIsNull && includeZero) {
                        partCountTest = std::pow(2.0, static_cast<double>(target - 1));
                    } else {
                        partCountTest = (includeZero) ? nChooseK(target + r - 1, r - 1) : nChooseK(target - 1, r - 1);
                    }
                }
                
                partCountTest = (bUserRows && numRows < partCountTest) ? numRows : partCountTest;
                
                if (partCountTest > std::numeric_limits<int>::max())
                    Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
                
                const std::size_t numParts = partCountTest;
                nCols = (isDistinct) ? distinctTest.limit : nCols;
                typeRcpp partitionsMatrix = (includeZero && (isComb || mIsNull)) ? 
                                                typeRcpp(numParts, nCols) : 
                                                    Rcpp::no_init_matrix(numParts, nCols);
                if (isDistinct) {
                    PartsStdDistinct(uR, numParts, partitionsMatrix, distinctTest.getAll,
                                     boundary, edge, lastCol, z, isComb, includeZero, mIsNull,  partitionsVec);
                } else {
                    PartitionsStandard(uR, numParts, partitionsMatrix, z, isComb,
                                       boundary, edge, lastCol, includeZero, mIsNull);
                }
                
                if (xtraCol)
                    for (std::size_t i = 0; i < numParts; ++i)
                        partitionsMatrix(i, r) = target;
                
                return partitionsMatrix;
            }
        }
        
        const double upperBound = std::min(vecMax, static_cast<double>(std::numeric_limits<int>::max()));
        const int maxRows = std::min(upperBound, numRows);
        
        if (bUserRows)
            partitionsVec.reserve(maxRows * r);
        
        int result = 0;
        
        if (isRep) {
            result = PartitionsRep(r, v, target, lastElem, lastCol, maxRows, isComb, partitionsVec);
        } else if (isMult) {
            result = PartitionsMultiSet(r, v, target, lastElem, lastCol, maxRows, isComb, Reps, partitionsVec);
        } else {
            result = PartitionsDistinct(r, v, target, lastElem, lastCol, maxRows, isComb, partitionsVec);
        }
        
        if (result) {
            std::size_t partitionLen = partitionsVec.size();
            std::size_t numResult = partitionLen / r;
            typeRcpp partitionsMatrix = Rcpp::no_init_matrix(numResult, nCols);
            
            for (std::size_t i = 0, k = 0; i < numResult; ++i)
                for (int j = 0; j < r; ++j, ++k)
                    partitionsMatrix(i, j) = partitionsVec[k];
                
            if (xtraCol)
                for (std::size_t i = 0; i < numResult; ++i)
                    partitionsMatrix(i, r) = target;
                    
            if (partitionLen >= upperBound) {
                Rcpp::warning("The algorithm terminated early as the number of results "
                              "meeting the criteria exceeds the container's maximum "
                              "capacity or 2^31 - 1");
            }
            
            return partitionsMatrix;
        } else {
            typeRcpp trivialRet(0, r);
            return trivialRet;
        }
    }
}

#endif
