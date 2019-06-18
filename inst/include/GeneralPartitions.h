#ifndef GENERAL_PARTITIONS_H
#define GENERAL_PARTITIONS_H

#include "CombPermUtils.h"

namespace Partitions {

    constexpr int solnExists = 1;
    constexpr int noSoln = 0;
    
    template <typename typeVector>
    int PartitionsRep(int n, int r, const std::vector<typeVector> &v, typeVector target,
                      double maxRows, bool isComb, double tol, std::vector<typeVector> &partitionsVec) {
    
        std::size_t count = 0;
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
    
        int numIter;
        
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
            
            int vertex = edge + 1;
            
            if (z[maxIndex] - z[edge] == 2) {
                vertex = maxIndex;
            } else {
                while (vertex < lastCol && (z[vertex] - z[edge]) < 2)
                    ++vertex;
            }
            
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
            
            if (vertex < maxIndex || z[maxIndex] < lastElem) {
                
                if (z[vertex] == z[edge])
                    ++vertex;
                
                while (vertex < pivot) {
                    int diVert = z[vertex] - z[edge];
                    int diPiv = lastElem - z[pivot];
                    
                    if (diVert == diPiv) {
                        z[vertex] -= diVert;
                        z[pivot] += diVert;
                        
                        ++vertex;
                        --pivot;
                    } else if (diVert < diPiv) {
                        z[vertex] -= diVert;
                        z[pivot] += diVert;
                        
                        ++vertex;
                    } else {
                        z[vertex] -= diPiv;
                        z[pivot] += diPiv;
                        
                        --pivot;
                    }
                }
                
                maxIndex = pivot;
                
                if (z[pivot] == lastElem) {
                    --pivot;
                } else if (pivot < lastCol && z[pivot] < z[pivot + 1]) {
                    ++maxIndex;
                }
            }
            
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
            
            if (vertex < outside || z[pivot] == outside - 1) {
                
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
                
                if (z[pivot] == lastElem) {
                    --pivot;
                } else if (pivot < lastCol && z[pivot + 1] - z[pivot] > 1) {
                    ++outside;
                }
            }
            
            while (outside > 0 && (z[outside] - z[outside - 1]) < 2)
                --outside;
            
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
        
        std::vector<typeVector> partitionsVec;
        const std::size_t calcRows = partitionsVec.max_size() / r;
        const std::size_t upperBound = std::min(calcRows, static_cast<std::size_t>(std::numeric_limits<int>::max()));
        const std::size_t maxRows = std::min(upperBound, static_cast<std::size_t>(numRows));
        
        if (bUserRows)
            partitionsVec.reserve(maxRows * r);
        
        std::sort(v.begin(), v.end());
        int result = noSoln;
        std::size_t nCols = (xtraCol) ? r + 1 : r;
        
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
