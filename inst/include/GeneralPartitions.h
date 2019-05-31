#ifndef GENERAL_PARTITIONS_H
#define GENERAL_PARTITIONS_H

#include "CombPermUtils.h"

template <typename typeVector>
struct valsData {
    std::vector<typeVector> v;
    std::vector<std::size_t> ind;
};

template <typename typeVector>
valsData<typeVector> valsExtreme(int n, int m, std::vector<typeVector> v, bool IsRep) {
    
    // When valsExtreme is called, we have already sorted v and myReps
    // in GeneralPartitions, otherwise we would need to do so below
    
    const unsigned long int uN = n;
    std::vector<typeVector> vals;
    std::vector<std::size_t> myIndex;
    
    if (IsRep) {
        const unsigned long int uM = m;
        vals.resize(v.size() * m);
        myIndex.resize(v.size() * m);
        
        for (std::size_t i = 0, k = 0; i < uN; ++i) {
            for (std::size_t j = 0; j < uM; ++j, ++k) {
                vals[k] = v[i];
                myIndex[k] = i;
            }
        }
        
    } else {
        vals.assign(v.cbegin(), v.cend());
        myIndex.resize(vals.size());
        std::iota(myIndex.begin(), myIndex.end(), 0);
    }
    
    valsData<typeVector> myVecs = {vals, myIndex};
    return myVecs;
}

template <typename typeVector>
std::vector<int> findStart(int n, int m, std::vector<typeVector> v,
                           bool IsRep, typeVector target, double tolerance) {
    
    std::vector<int> valInd(m);
    const unsigned long int lastElem = m - 1;
    const typeVector zero = static_cast<typeVector>(0);
    
    valsData<typeVector> myVecs = valsExtreme(n, m, v, IsRep);
    std::vector<typeVector> vals = myVecs.v;
    std::vector<std::size_t> myIndex = myVecs.ind;
    const unsigned long int valSize = vals.size();
    
    // if the largest sum is smaller than target, there is no solution
    typeVector testMax = std::accumulate(vals.cend() - m, vals.cend(), zero);
    const std::vector<int> noSoln(1, 0);
    
    if (testMax <= target)  {
        if (testMax == target) {
            std::vector<int> maxSoln;
            
            for (auto it = myIndex.cend() - m; it < myIndex.cend(); ++it)
                maxSoln.push_back(*it);
            
            return maxSoln;
        }
        
        return noSoln;
    }
    
    // if the smallest sum is greater than target, there is no solution
    typeVector testMin = std::accumulate(vals.cbegin(), vals.cbegin() + m, zero);
    
    if (testMin >= target)  {
        if (testMin == target) {
            std::vector<int> minSoln;
            
            for (auto it = myIndex.cbegin(); it < myIndex.cbegin() + m; ++it)
                minSoln.push_back(*it);
            
            return minSoln;
        }
        
        return noSoln;
    }
    
    int nextPos = m;
    int mid = (valSize - nextPos) / 2;
    
    // Perform a modified binary search. N.B. Since vals is sorted
    // in ascending order, we must have the last m - 1 elements.
    typeVector partial = std::accumulate(vals.cend() - (m - 1), vals.cend(), zero);
    typeVector dist = target - (partial + vals[mid]);
    
    int lowBnd = (dist > 0) ? mid : 0;
    int uppBnd = (dist > 0) ? valSize - nextPos : mid;
    int ind = 0;
    std::size_t i = 0, prevInd = 0;
    
    for (; i <= lastElem; ++i) {
        while ((uppBnd - lowBnd) > 1 && dist != 0) {
            mid = (uppBnd - lowBnd) / 2;
            ind = lowBnd + mid;
            dist = target - (partial + vals[ind]);
            
            if (dist > 0)
                lowBnd = ind;
            else
                uppBnd = ind;
        }
        
        if (ind > 0)
            while (vals[ind] == vals[ind - 1] && ind > prevInd)
                --ind;
        
        // Check last index. N.B. There are some cases when
        // ind == lowBnd and dist < 0. This will not matter
        // as we simply reassign ind and recompute dist
        if (dist < 0) {
            ind = lowBnd;
            dist = target - (partial + vals[ind]);
        }
        
        // We must have dist < 0. Below is an informal proof.
        // The sub-sequences are defined as below:
        //     A_max = a_(i + 1), a_(i + 2), ..., a_m
        //     A_set = a_1, a_2, ..., a_(i - 1)
        // A_set are those elements that have already been
        // determined by the algorithm. A_max is maximal (i.e. 
        // constructed of the the (m - i) largest elements).
        // We seek to determine the i_th element given the
        // following contraints:
        //         A_sum = A_set + a_i + A_max
        //         dist = target - A_sum
        // With the goal of finding the minimum lexicographic
        // combination such that the dist = 0 (i.e. target =
        // A_sum). If we have dist > 0 for any i, then it will
        // be impossible to obtain dist = 0. dist > 0 implies
        // that the target > A_sum, and since A_max is already
        // maximal, we are not able to increase A_sum in later
        // iterations, thus we must have dist <= 0 for all i.
        if (dist > 0 && ind < valSize)
            ++ind;
        
        // We must ensure that we have the smallest index, ind,
        // of vals such that ind > valInd[i - 1] and vals[ind]
        // remains unchanged
        
        // We found our index
        prevInd = valInd[i] = ind;
        partial += vals[ind];
        --nextPos;
        
        if (ind < valSize)
            ++ind;
        
        // We adjust the length that we need to span as well as
        // reset our lower and upper bounds.
        mid = lowBnd = ind;
        uppBnd = valSize - nextPos;
        
        if (i < lastElem) {
            partial -= vals[uppBnd];
            dist = target - (partial + vals[mid]);
        }
    }
    
    std::vector<int> res(m);
    std::vector<typeVector> check(m);
    
    for (std::size_t i = 0; i <= lastElem; ++i) {
        res[i] = myIndex[valInd[i]];
        check[i] = vals[valInd[i]];
    }
    
    // The algorithm above finds the first possible sum that equals
    // target. If there is no combination of elements from v that sum
    // to target, the algo returns the combination such that its sum
    // is closest to target and greater than target
    typeVector finalCheck = std::accumulate(check.cbegin(), check.cend(), zero);
    
    if (!std::is_integral<typeVector>::value) {
        if (std::abs(finalCheck - target) > tolerance)
            return noSoln;
    } else if (finalCheck != target) {
        return noSoln;
    }
    
    return res;
}

template <typename typeRcpp, typename typeVector>
typeRcpp GeneralPartitions(int n, int r, std::vector<typeVector> &v, bool isRep, typeVector target,
                           double numRows, bool isComb, bool xtraCol, bool bUserRows, double tol) {
    
    std::size_t count = 0;
    std::vector<typeVector> partitionsVec;
    const std::size_t maxRows = std::min(static_cast<double>(
        std::numeric_limits<int>::max()), numRows);
    
    std::sort(v.begin(), v.end());
    std::vector<int> z = findStart(n, r, v, isRep, target, tol);
    
    // Check to ensure there exist at least one solution
    if (z.size() == 1) {
        if (z[0] == 0) {
            typeRcpp noSoln;
            return noSoln;
        }
    }
    
    std::vector<typeVector> testVec(r);
    std::vector<int> zCheck;
    int numIter, maxZ = n - 1;
    const int r1 = r - 1;
    
    if (bUserRows)
        partitionsVec.reserve(maxRows * r);
    
    if (isRep) {
        // smallest index such that z[maxIndex] == currMax
        int maxIndex = r1;
        int currMax = z[maxIndex];
        
        while (maxIndex > 0 && z[maxIndex - 1] == currMax)
            --maxIndex;
        
        // pivot is the greatest index such that z[pivot] < maxZ
        // We know that if z[maxIndex] < maxZ ==>> pivot = r1
        int pivot = (z[maxIndex] == maxZ) ? maxIndex - 1 : r1;
        
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
                while (vertex < r1 && (z[vertex] - z[edge]) < 2)
                    ++vertex;
            }
            
            ++z[edge];
            --z[vertex];
            
            if (vertex == maxIndex) {
                if (maxIndex < r1)
                    ++maxIndex;
                
                currMax = z[maxIndex];
                
                while (maxIndex > 0 && z[maxIndex - 1] == currMax)
                    --maxIndex;
                
                pivot = (z[maxIndex] == maxZ) ? maxIndex - 1 : r1;
            }
            
            if (vertex < maxIndex || z[maxIndex] < maxZ) {
                
                if (z[vertex] == z[edge])
                    ++vertex;
                
                while (vertex < pivot) {
                    int diVert = z[vertex] - z[edge];
                    int diPiv = maxZ - z[pivot];
                    
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
                
                if (z[pivot] == maxZ) {
                    --pivot;
                } else if (pivot < r1 && z[pivot] < z[pivot + 1]) {
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
        }
        
    } else {
        
        int indexRows = isComb ? 0 : static_cast<int>(NumPermsNoRep(r, r1));
        auto indexMatrix = std::make_unique<int[]>(indexRows * r);
        
        if (!isComb) {
            indexRows = static_cast<int>(NumPermsNoRep(r, r1));
            std::vector<int> indexVec(r);
            std::iota(indexVec.begin(), indexVec.end(), 0);
            
            for (int i = 0, myRow = 0; i < indexRows; ++i, myRow += r) {
                for (int j = 0; j < r; ++j)
                    indexMatrix[myRow + j] = indexVec[j];
                
                std::next_permutation(indexVec.begin(), indexVec.end());
            }
        }
        
        // Largest index such that z[outside] - z[outside - 1] > 1
        int outside = r1;
        
        while (outside > 0 && (z[outside] - z[outside - 1]) < 2)
            --outside;
        
        // pivot is the greatest index that can be incremented...
        // Either z[pivot + 1] - z[pivot] > 1 or if z[r1] < maxZ
        // pivot = r1 since incrementing z[r1] is possible
        int pivot = (z[r1] < maxZ) ? r1 : outside - 1;
        
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
                        partitionsVec.push_back(v[z[k]]);
                
                count += indexRows;
            }
            
            if (count >= maxRows)
                break;
            
            int vertex = edge + 1;
            tarDiff = 3;
            
            while (vertex < r1 && (z[vertex] - z[edge]) < tarDiff) {
                ++vertex;
                ++tarDiff;
            }
            
            ++z[edge];
            --z[vertex];
            
            if (vertex == outside) {
                if (outside < r1)
                    ++outside;
                
                while (outside > 0 && (z[outside] - z[outside - 1]) < 2)
                    --outside;
                
                pivot = (z[r1] < maxZ) ? r1 : outside - 1;
            }
            
            if (vertex < outside || z[pivot] == outside - 1) {
                
                if (z[vertex] - z[vertex - 1] == 1)
                    ++vertex;
                
                while (vertex < pivot) {
                    --z[vertex];
                    ++z[pivot];
                    
                    if (z[vertex] - z[vertex - 1] == 1)
                        ++vertex;
                    
                    if ((pivot < r1 && z[pivot + 1] - z[pivot] == 1) || z[pivot] == maxZ)
                        --pivot;
                }
                
                outside = pivot;
                
                if (z[pivot] == maxZ) {
                    --pivot;
                } else if (pivot < r1 && z[pivot + 1] - z[pivot] > 1) {
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
                
                ++count;
            } else {
                if (indexRows + count > maxRows)
                    indexRows = maxRows - count;
                
                for (int j = 0, myRow = 0; j < indexRows; ++j, myRow += r)
                    for (int k = 0; k < r; ++k)
                        partitionsVec.push_back(v[z[k]]);
                
                count += indexRows;
            }
        }
    }
    
    unsigned long int nCols = (xtraCol) ? r + 1 : r;
    typeRcpp partitionsMatrix = Rcpp::no_init_matrix(count, nCols);
    
    for (std::size_t i = 0, k = 0; i < count; ++i)
        for (int j = 0; j < r; ++j, ++k)
            partitionsMatrix(i, j) = partitionsVec[k];
    
    if (xtraCol)
        for (std::size_t i = 0; i < count; ++i)
            partitionsMatrix(i, r) = target;
    
    if (count > std::numeric_limits<int>::max()) {
        Rcpp::warning("The algorithm terminated early as the number of "
                          "results meeting the criteria exceeds 2^31 - 1.");
    }
    
    return partitionsMatrix;
}

#endif
