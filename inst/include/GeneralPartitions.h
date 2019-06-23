#ifndef GENERAL_PARTITIONS_H
#define GENERAL_PARTITIONS_H

#include "CombPermUtils.h"

namespace Partitions {

    constexpr int solnExists = 1;
    constexpr int noSoln = 0;
    
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
                for (int j = 0; j < n1 - i + 1; ++j)
                    for (int k = 0; k < i; ++k)
                        memoize(j, k) = 0;

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
    
    template <typename typeRcpp>
    void PartitionsStandard(int n, int r, int numRows, bool isComb, 
                            bool includeZero, typeRcpp &partitionsMatrix) {
        
        const int lastElem = n - 1;
        const int lastCol = r - 1;
        
        std::vector<int> z(r, 0);
        z[lastCol] = lastElem;
        
        if (!includeZero) {
            std::fill(z.begin(), z.end(), 1);
            z[lastCol] = n - r + 1;
        }
        
        int numIter = 0;
        int count = 0;
        int limitRows = numRows - 1;
        int maxIndex = lastCol;
        int pivot = (includeZero) ? lastCol - 1 : lastCol;
        int edge = maxIndex - 1;
        
        for (; count < limitRows;) {
            if (isComb) {
                for (int k = 0; k < r; ++k)
                    partitionsMatrix(count, k) = z[k];
                
                ++count;
            } else {
                numIter = static_cast<int>(NumPermsWithRep(z));
                
                if ((numIter + count) > numRows)
                    numIter = numRows - count;
                
                for (int i = 0; i < numIter; ++i, ++count) {
                    for (int k = 0; k < r; ++k)
                        partitionsMatrix(count, k) = z[k];

                    std::next_permutation(z.begin(), z.end());
                }
            }
            
            int vertex = (z[maxIndex] - z[edge] == 2) ? maxIndex : edge + 1;
            
            ++z[edge];
            --z[vertex];
            
            if (vertex == maxIndex) {
                if (maxIndex < lastCol)
                    ++maxIndex;
                
                int currMax = z[maxIndex];
                
                while (z[maxIndex - 1] == currMax)
                    --maxIndex;
                
                pivot = (z[maxIndex] == lastElem) ? maxIndex - 1 : lastCol;
            }
            
            if (z[vertex] == z[edge])
                ++vertex;
            
            while (vertex < pivot) {
                int distVert = z[vertex] - z[edge];
                int distPivot = lastElem - z[pivot];
                
                if (distVert == distPivot) {
                    z[vertex] -= distVert;
                    z[pivot] += distVert;
                    
                    ++vertex;
                    --pivot;
                } else if (distVert < distPivot) {
                    z[vertex] -= distVert;
                    z[pivot] += distVert;
                    
                    ++vertex;
                } else {
                    z[vertex] -= distPivot;
                    z[pivot] += distPivot;
                    
                    --pivot;
                }
            }
            
            maxIndex = pivot;
            
            if (pivot < lastCol && z[pivot] < z[pivot + 1])
                ++maxIndex;
            
            int currMax = z[maxIndex];
            
            while (z[maxIndex - 1] == currMax)
                --maxIndex;
            
            edge = maxIndex - 1;
            
            while (z[maxIndex] - z[edge] < 2)
                --edge;
        }
        
        if (isComb) {
            for (int k = 0; k < r; ++k)
                partitionsMatrix(limitRows, k) = z[k];
        } else {
            numIter = static_cast<int>(NumPermsWithRep(z));
            
            if ((numIter + count) > numRows)
                numIter = numRows - count;
            
            for (int i = 0; i < numIter; ++i, ++count) {
                for (int k = 0; k < r; ++k)
                    partitionsMatrix(count, k) = z[k];

                std::next_permutation(z.begin(), z.end());
            }
        }
    }
    
    template <typename typeVector>
    int PartitionsRep(int n, int r, const std::vector<typeVector> &v, typeVector target,
                      double maxRows, bool isComb, double tol, std::vector<typeVector> &partitionsVec) {
    
        std::vector<int> z(r);
        const int lastElem = n - 1;
        const int lastCol = r - 1;
        
        typeVector testMax = v[lastElem] * r;
        if (testMax < target)  {return noSoln;}
        
        typeVector partial = testMax;
        partial -= v[lastElem];
        
        typeVector testMin = v[0] * r;
        if (testMin > target)  {return noSoln;}
        
        int mid = lastElem / 2;
        typeVector dist = target - (partial + v[mid]);
        
        int lowBnd = (dist > 0) ? mid : 0;
        int uppBnd = (dist > 0) ? lastElem : mid;
        int ind = mid;
        
        for (int i = 0; i < r; ++i) {
            while ((uppBnd - lowBnd) > 1 && dist != 0) {
                mid = (uppBnd - lowBnd) / 2;
                ind = lowBnd + mid;
                dist = target - (partial + v[ind]);
                
                if (dist > 0)
                    lowBnd = ind;
                else
                    uppBnd = ind;
            }
            
            // Check last index. N.B. There are some cases when
            // ind == lowBnd and dist < 0. This will not matter
            // as we simply reassign ind and recompute dist
            if (dist < 0) {
                ind = lowBnd;
                dist = target - (partial + v[ind]);
            }
            
            if (dist > 0 && ind < lastElem)
                ++ind;
            
            z[i] = ind;
            partial += v[ind];
            
            lowBnd = ind;
            uppBnd = lastElem;
            mid = (uppBnd - lowBnd) / 2;
            
            ind = lowBnd + mid;
            partial -= v[lastElem];
            dist = target - (partial + v[ind]);
        }
        
        std::vector<typeVector> check(r);
        
        for (int i = 0; i < r; ++i)
            check[i] = v[z[i]];
        
        // The algorithm above finds the first possible sum that equals
        // target. If there is no combination of elements from v that sum
        // to target, the algo returns the combination such that its sum
        // is closest to target and greater than target
        typeVector finalCheck = std::accumulate(check.cbegin(), check.cend(), static_cast<typeVector>(0));
    
        if (!std::is_integral<typeVector>::value) {
            if (std::abs(finalCheck - target) > tol)
                return noSoln;
        } else if (finalCheck != target) {
            return noSoln;
        }
    
        int numIter = 0;
        std::size_t count = 0;
        
        // smallest index such that z[maxIndex] == currMax
        int maxIndex = lastCol;
        int currMax = z[maxIndex];
        
        while (maxIndex > 0 && z[maxIndex - 1] == currMax)
            --maxIndex;
        
        // pivot is the greatest index such that z[pivot] < lastElem
        // We know that if z[maxIndex] < lastElem ==>> pivot = lastCol
        int pivot = (z[maxIndex] == lastElem) ? maxIndex - 1 : lastCol;
        
        // edge is the greatest index such that z[maxIndex] - z[edge] >= 2
        // This is the index that will be be used as a starting point
        // to determine the next combination that meets the criteria
        int edge = maxIndex - 1;
        
        while (edge > 0 && (z[maxIndex] - z[edge]) < 2)
            --edge;
        
        while (edge >= 0 && (z[maxIndex] - z[edge]) >= 2) {
            if (isComb) {
                for (int k = 0; k < r; ++k)
                    partitionsVec.push_back(v[z[k]]);
                
                ++count;
            } else {
                numIter = static_cast<int>(NumPermsWithRep(z));
                
                if ((numIter + count) > maxRows)
                    numIter = maxRows - count;
                
                for (int i = 0; i < numIter; ++i) {
                    for (int k = 0; k < r; ++k)
                        partitionsVec.push_back(v[z[k]]);
                    
                    std::next_permutation(z.begin(), z.end());
                }
                
                count += numIter;
            }
            
            if (count >= maxRows)
                break;
            
            int vertex = (z[maxIndex] - z[edge] == 2) ? maxIndex : edge + 1;
            
            ++z[edge];
            --z[vertex];
            
            if (vertex == maxIndex) {
                if (maxIndex < lastCol)
                    ++maxIndex;
                
                currMax = z[maxIndex];
                
                while (maxIndex > 0 && z[maxIndex - 1] == currMax)
                    --maxIndex;
                
                pivot = (z[maxIndex] == lastElem) ? maxIndex - 1 : lastCol;
            }
            
            if (z[vertex] == z[edge])
                ++vertex;
            
            while (vertex < pivot) {
                int distVert = z[vertex] - z[edge];
                int distPivot = lastElem - z[pivot];
                
                if (distVert == distPivot) {
                    z[vertex] -= distVert;
                    z[pivot] += distVert;
                    
                    ++vertex;
                    --pivot;
                } else if (distVert < distPivot) {
                    z[vertex] -= distVert;
                    z[pivot] += distVert;
                    
                    ++vertex;
                } else {
                    z[vertex] -= distPivot;
                    z[pivot] += distPivot;
                    
                    --pivot;
                }
            }
                
            maxIndex = pivot;
            
            if (pivot < lastCol && z[pivot] < z[pivot + 1])
                ++maxIndex;
            
            currMax = z[maxIndex];
            
            while (maxIndex > 0 && z[maxIndex - 1] == currMax)
                --maxIndex;
            
            edge = maxIndex - 1;
            
            while (edge > 0 && (z[maxIndex] - z[edge]) < 2)
                --edge;
        }
        
        if (count < maxRows) {
            if (isComb) {
                for (int k = 0; k < r; ++k)
                    partitionsVec.push_back(v[z[k]]);
            } else {
                numIter = static_cast<int>(NumPermsWithRep(z));
                
                if ((numIter + count) > maxRows)
                    numIter = maxRows - count;
                
                for (int i = 0; i < numIter; ++i) {
                    for (int k = 0; k < r; ++k)
                        partitionsVec.push_back(v[z[k]]);
                    
                    std::next_permutation(z.begin(), z.end());
                }
            }
        }
        
        return solnExists;
    }
    
    template <typename typeVector>
    int PartitionsDistinct(int n, int r, const std::vector<typeVector> &v, typeVector target,
                           double maxRows, bool isComb, double tol, std::vector<typeVector> &partitionsVec) {
        
        std::size_t count = 0;
        std::vector<int> z(r);
        const int lastElem = n - 1;
        const int lastCol = r - 1;
        
        typeVector testMax = std::accumulate(v.cend() - r, v.cend(), static_cast<typeVector>(0));
        if (testMax < target)  {return noSoln;}
        
        int currPos = n - r;
        typeVector partial = testMax;
        partial -= v[currPos];
        
        typeVector testMin = std::accumulate(v.cbegin(), v.cbegin() + r, static_cast<typeVector>(0));
        if (testMin > target)  {return noSoln;}
        
        int mid = currPos / 2;
        int dist = target - (partial + v[mid]);
        
        int lowBnd = (dist > 0) ? mid : 0;
        int uppBnd = (dist > 0) ? currPos : mid;
        int ind = mid;
        
        for (int i = 0; i < r; ++i) {
            while ((uppBnd - lowBnd) > 1 && dist != 0) {
                mid = (uppBnd - lowBnd) / 2;
                ind = lowBnd + mid;
                dist = target - (partial + v[ind]);
                
                if (dist > 0)
                    lowBnd = ind;
                else
                    uppBnd = ind;
            }
            
            // Check last index. N.B. There are some cases when
            // ind == lowBnd and dist < 0. This will not matter
            // as we simply reassign ind and recompute dist
            if (dist < 0) {
                ind = lowBnd;
                dist = target - (partial + v[ind]);
            }
            
            if (dist > 0 && ind < lastElem)
                ++ind;
            
            z[i] = ind;
            partial += v[ind];
            
            ++ind;
            ++currPos;
            
            lowBnd = ind;
            uppBnd = currPos;
            mid = (uppBnd - lowBnd) / 2;
            
            ind = lowBnd + mid;
            partial -= v[currPos];
            dist = target - (partial + v[ind]);
        }
        
        std::vector<typeVector> check(r);
        
        for (int i = 0; i < r; ++i)
            check[i] = v[z[i]];
        
        // The algorithm above finds the first possible sum that equals
        // target. If there is no combination of elements from v that sum
        // to target, the algo returns the combination such that its sum
        // is closest to target and greater than target
        typeVector finalCheck = std::accumulate(check.cbegin(), check.cend(), static_cast<typeVector>(0));
        
        if (!std::is_integral<typeVector>::value) {
            if (std::abs(finalCheck - target) > tol)
                return noSoln;
        } else if (finalCheck != target) {
            return noSoln;
        }
        
        int indexRows = isComb ? 0 : static_cast<int>(NumPermsNoRep(r, lastCol));
        auto indexMatrix = std::make_unique<int[]>(indexRows * r);
        
        if (!isComb) {
            indexRows = static_cast<int>(NumPermsNoRep(r, lastCol));
            std::vector<int> indexVec(r);
            std::iota(indexVec.begin(), indexVec.end(), 0);
            
            for (int i = 0, myRow = 0; i < indexRows; ++i, myRow += r) {
                for (int j = 0; j < r; ++j)
                    indexMatrix[myRow + j] = indexVec[j];
                
                std::next_permutation(indexVec.begin(), indexVec.end());
            }
        }
        
        // Largest index such that z[outside] - z[outside - 1] > 1
        int outside = lastCol;
        
        while (outside > 0 && (z[outside] - z[outside - 1]) < 2)
            --outside;
        
        // pivot is the greatest index that can be incremented...
        // Either z[pivot + 1] - z[pivot] > 1 or if z[lastCol] < lastElem
        // pivot = lastCol since incrementing z[lastCol] is possible
        int pivot = (z[lastCol] < lastElem) ? lastCol : outside - 1;
        
        // edge is the greatest index such that when incremented
        // the result will be at least one less than its neighbor
        // even if its neighbor is decremented
        int edge = outside - 1;
        int tarDiff = 3;
        
        while (edge > 0 && (z[outside] - z[edge]) < tarDiff) {
            --edge;
            ++tarDiff;
        }

        while (edge >= 0 && (z[outside] - z[edge]) >= tarDiff) {
            if (isComb) {
                for (int k = 0; k < r; ++k)
                    partitionsVec.push_back(v[z[k]]);
                
                ++count;
            } else {
                if (indexRows + count > maxRows)
                    indexRows = maxRows - count;
                
                for (int j = 0, myRow = 0; j < indexRows; ++j, myRow += r)
                    for (int k = 0; k < r; ++k)
                        partitionsVec.push_back(v[z[indexMatrix[myRow + k]]]);
                
                count += indexRows;
            }
            
            if (count >= maxRows)
                break;
            
            int vertex = edge + 1;
            tarDiff = 3;
            
            while (vertex < lastCol && (z[vertex] - z[edge]) < tarDiff) {
                ++vertex;
                ++tarDiff;
            }
            
            ++z[edge];
            --z[vertex];
            
            if (vertex == outside) {
                if (outside < lastCol)
                    ++outside;
                
                while (outside > 0 && (z[outside] - z[outside - 1]) < 2)
                    --outside;
                
                pivot = (z[lastCol] < lastElem) ? lastCol : outside - 1;
            }
            
            if (vertex < outside) {
                if (z[vertex] - z[vertex - 1] == 1)
                    ++vertex;
                
                while (vertex < pivot) {
                    --z[vertex];
                    ++z[pivot];
                    
                    if (z[vertex] - z[vertex - 1] == 1)
                        ++vertex;
                    
                    if ((pivot < lastCol && z[pivot + 1] - z[pivot] == 1) || z[pivot] == lastElem)
                        --pivot;
                }
                
                outside = pivot;
                
                if (pivot < lastCol && z[pivot + 1] - z[pivot] > 1)
                    ++outside;
            }
            
            edge = outside - 1;
            tarDiff = 3;
            
            while (edge > 0 && (z[outside] - z[edge]) < tarDiff) {
                --edge;
                ++tarDiff;
            }
        }
        
        if (count < maxRows) {
            if (isComb) {
                for (int k = 0; k < r; ++k)
                    partitionsVec.push_back(v[z[k]]);
            } else {
                if (indexRows + count > maxRows)
                    indexRows = maxRows - count;
                
                for (int j = 0, myRow = 0; j < indexRows; ++j, myRow += r)
                    for (int k = 0; k < r; ++k)
                        partitionsVec.push_back(v[z[indexMatrix[myRow + k]]]);
            }
        }
        
        return solnExists;
    }
    
    template <typename typeRcpp, typename typeVector>
    typeRcpp GeneralPartitions(int n, int r, std::vector<typeVector> &v, typeVector target, bool isRep,
                               double numRows, bool isComb, bool xtraCol, bool bUserRows, double tol) {
        
        std::sort(v.begin(), v.end());
        typeVector myMax = v.back();
        std::size_t nCols = (xtraCol) ? r + 1 : r;
        
        if (myMax == target && (v.size() == target || v.size() == (target + 1)) && isRep && r == target) {
            const bool includeZero = (v.front() == 1) ? false : true;
            const double countTest = partitionCount(n, r, includeZero);
            
            if (countTest > std::numeric_limits<int>::max())
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            const std::size_t numParts = (bUserRows && numRows < countTest) ? numRows : countTest;
            typeRcpp partitionsMatrix = Rcpp::no_init_matrix(numParts, nCols);
            PartitionsStandard(n, r, numParts, isComb, includeZero, partitionsMatrix);
            
            if (xtraCol)
                for (std::size_t i = 0; i < numParts; ++i)
                    partitionsMatrix(i, r) = target;
            
            return partitionsMatrix;
        }
        
        std::vector<typeVector> partitionsVec;
        const double calcRows = partitionsVec.max_size() / r;
        const double upperBound = std::min(calcRows, static_cast<double>(std::numeric_limits<int>::max()));
        const double maxRows = std::min(static_cast<double>(upperBound), numRows);
        
        if (bUserRows)
            partitionsVec.reserve(maxRows * r);
        
        int result = noSoln;
        
        if (isRep)
            result = PartitionsRep(n, r, v, target, maxRows, isComb, tol, partitionsVec);
        else
            result = PartitionsDistinct(n, r, v, target, maxRows, isComb, tol, partitionsVec);
        
        if (result) {
            std::size_t partitionLen = partitionsVec.size();
            
            if ((partitionLen % r) != 0)
                Rcpp::stop("Unexpected number of results");
            
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
