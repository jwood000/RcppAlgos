#include "CombPermUtils.h"
#include "Cpp14MakeUnique.h"


void BinaryNextElem(int &uppBnd, int &lowBnd, int &ind, int lastElem,
                           int64_t target, int64_t partial, const std::vector<int64_t> &v) {
    
    int64_t dist = target - (partial + v[ind]);
    
    while ((uppBnd - lowBnd) > 1 && dist != 0) {
        const int mid = (uppBnd - lowBnd) / 2;
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
    
    if (dist > 0 && ind < lastElem)
        ++ind;
}

int GetFirstCombo(int r, const std::vector<int64_t> &v, bool isRep, bool isMult,
                  std::vector<int> &z, const std::vector<int> &freqs, int64_t target,
                  const std::vector<int> &Reps, int lastCol, int lastElem) {
    
    int64_t testMax = 0;
    constexpr int64_t zero64 = 0;
    
    if (isRep) {
        testMax = v[lastElem] * r;
    } else if (isMult) {
        const int lenMinusR = freqs.size() - r;
        
        for (int i = freqs.size() - 1, j = 0; i >= lenMinusR; --i, ++j)
            testMax += v[freqs[i]];
    } else {
        testMax = std::accumulate(v.cend() - r, v.cend(), zero64);
    }
    
    if (testMax < target)  {return 0;}
    int zExpCurrPos = (isMult) ? freqs.size() - r : 0;
    int currPos = (isMult) ? freqs[zExpCurrPos] : ((isRep) ? lastElem : (v.size() - r));
    
    int64_t partial = testMax;
    partial -= v[currPos];
    int64_t testMin = 0;
    
    if (isRep) {
        testMin = v[0] * r;
    } else if (isMult) {
        for (int i = 0; i < r; ++i)
            testMin += v[freqs[i]];
    } else {
        testMin = std::accumulate(v.cbegin(), v.cbegin() + r, zero64);
    }
    
    if (testMin > target)  {return 0;}
    int mid = currPos / 2;
    int64_t dist = target - (partial + v[mid]);
    
    int lowBnd = (dist > 0) ? mid : 0;
    int uppBnd = (dist > 0) ? currPos : mid;
    int ind = mid;
    
    std::vector<int> repsCounter;
    
    if (isMult)
        repsCounter.assign(Reps.cbegin(), Reps.cend());
    
    for (int i = 0; i < lastCol; ++i) {
        BinaryNextElem(uppBnd, lowBnd, ind, lastElem, target, partial, v);
        z[i] = ind;
        partial += v[ind];
        
        if (isMult) {
            --repsCounter[ind];
            
            if (repsCounter[ind] == 0)
                ++ind;
            
            ++zExpCurrPos;
            currPos = freqs[zExpCurrPos];
        } else if (!isRep) {
            ++ind;
            ++currPos;
        }
        
        lowBnd = ind;
        uppBnd = currPos;
        mid = (uppBnd - lowBnd) / 2;
        
        ind = lowBnd + mid;
        partial -= v[currPos];
    }
    
    BinaryNextElem(uppBnd, lowBnd, ind, lastElem, target, partial, v);
    z[lastCol] = ind;
    
    // The algorithm above finds the first possible sum that equals
    // target. If there is no combination of elements from v that sum
    // to target, the algo returns the combination such that its sum
    // is closest to target and greater than target
    int64_t finalCheck = 0;
    
    for (int i = 0; i < r; ++i)
        finalCheck += v[z[i]];
    
    if (finalCheck != target)
        return 0;
    
    return 1;
}

inline void PopulateVec(int r, const std::vector<int64_t> &v, std::vector<int> &z,
                        int &count, int maxRows, bool isComb, std::vector<int64_t> &partitionsVec) {
    if (count < maxRows) {
        if (isComb) {
            for (int k = 0; k < r; ++k)
                partitionsVec.push_back(v[z[k]]);
            
            ++count;
        } else {
            do {
                for (int k = 0; k < r; ++k)
                    partitionsVec.push_back(v[z[k]]);
                
                ++count;
            } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
        }
    }
}

inline bool BndDecrementPossible(const std::vector<int> &rpsCnt,
                                 const std::vector<int> &z, int boundary) {
    if (boundary > 1) {
        if ((z[boundary] - z[boundary - 1]) < 2) {
            return ((z[boundary] != z[boundary - 1]) ? !rpsCnt[z[boundary] - 1] : true);
        } else {
            return false;
        }
    } else {
        return false;
    }
}

inline bool VtxDecrementPossible(const std::vector<int> &rpsCnt, int lastCol,
                                 const std::vector<int> &z, int vertex, int edge) {
    if (vertex < lastCol) {
        const int myDiff = z[vertex] - z[edge];
        
        if (myDiff < 2) {
            return true;
        } else  {
            return ((myDiff != 2 || rpsCnt[z[edge] + 1] <= 1) &&
                    (myDiff <= 2 || !rpsCnt[z[edge] + 1] || !rpsCnt[z[vertex] - 1]));
        }
    } else {
        return false;
    }
}

inline bool EdgeIncrementPossible(const std::vector<int> &rpsCnt,
                                  const std::vector<int> &z, int edge, int boundary) {
    if (edge) {
        const int myDiff = z[boundary] - z[edge];
        
        if (myDiff < 2) {
            return true;
        } else {
            return (myDiff != 2 || rpsCnt[z[edge] + 1] < 2) && (myDiff <= 2 || !rpsCnt[z[edge] + 1]);
        }
    } else {
        return false;
    }
}

inline int GetPivotExtr(const std::vector<int> &rpsCnt,
                        const std::vector<int> &z, int lastCol, int lastElem) {
    int res = lastCol - 1;
    
    while (res > 0 && z[res] == lastElem)
        --res;
    
    while (res > 0 && !rpsCnt[z[res] + 1])
        --res;
    
    return res;
}

inline bool PivotDecrementPossible(const std::vector<int> &rpsCnt, int lastElem,
                                   const std::vector<int> &z, int pivot, int vertex) {
    if (pivot > vertex) {
        if (z[pivot] == lastElem) {
            return true;
        } else {
            return !rpsCnt[z[pivot] + 1];
        }
    } else {
        return false;
    }
}

inline bool keepGoing(const std::vector<int> &rpsCnt, int lastElem,
                      const std::vector<int> &z, int edge, int boundary) {
    if (edge >= 0) {
        const int myDiff = z[boundary] - z[edge];
        
        if (myDiff < 2) {
            return false;
        } else if (myDiff == 2) {
            return (rpsCnt[z[edge] + 1] > 1);
        } else {
            return (rpsCnt[z[edge] + 1] && rpsCnt[z[boundary] - 1]);
        }
    } else {
        return false;
    }
}

int PartitionsMultiSet(int r, const std::vector<int64_t> &v, int64_t target, 
                       int lastElem, int lastCol, int maxRows, bool isComb,
                       const std::vector<int> &Reps, std::vector<int64_t> &partitionsVec) {
    
    std::vector<int> z(r);
    std::vector<int> zExpanded;
    
    for (std::size_t i = 0; i < v.size(); ++i)
        for (int j = 0; j < Reps[i]; ++j)
            zExpanded.push_back(i);
    
    if (!GetFirstCombo(r, v, false, true, z, zExpanded, target, Reps, lastCol, lastElem))
        return 0;
    
    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());
    
    for (const auto ind: z)
        --rpsCnt[ind];
    
    // boundary is the greatest index that can be decremented
    int b = lastCol;
    
    while (BndDecrementPossible(rpsCnt, z, b))
        --b;
    
    // pivot is the greatest index that can be incremented
    int p = (z[lastCol] < lastElem) ? lastCol : GetPivotExtr(rpsCnt, z, lastCol, lastElem);
    
    // edge is the greatest index smaller than boundary that can be incremented
    int e = b - 1;
    
    while (EdgeIncrementPossible(rpsCnt, z, e, b))
        --e;
    
    int count = 0;
    
    while (keepGoing(rpsCnt, lastElem, z, e, b)) {
        
        PopulateVec(r, v, z, count, maxRows, isComb, partitionsVec);
        
        if (count >= maxRows)
            break;
        
        // vertex is the smallest index greater than edge that can be decremented
        int v = e + 1;
        
        while (VtxDecrementPossible(rpsCnt, lastCol, z, v, e))
            ++v;
        
        ++rpsCnt[z[e]];
        ++z[e];
        --rpsCnt[z[e]];
        
        ++rpsCnt[z[v]];
        --z[v];
        --rpsCnt[z[v]];
        
        if (v == b) {
            if (b < lastCol)
                ++b;
            
            while (BndDecrementPossible(rpsCnt, z, b))
                --b;
            
            p = (z[lastCol] < lastElem) ? lastCol : GetPivotExtr(rpsCnt, z, lastCol, lastElem);
        }
        
        while ((v < lastCol) && (z[v] == z[v - 1] || 
               z[v] == z[e] || (z[v] - z[v - 1] == 1 && !rpsCnt[z[v - 1]]))) {
            ++v;
        }
        
        while (v < p && rpsCnt[z[v] - 1] && rpsCnt[z[p] + 1]) {
            ++rpsCnt[z[v]];
            --z[v];
            --rpsCnt[z[v]];
            
            ++rpsCnt[z[p]];
            ++z[p];
            --rpsCnt[z[p]];
            
            while (z[v] == z[v - 1] || (z[v] - z[v - 1] == 1 && !rpsCnt[z[v - 1]]))
                ++v;
            
            while (PivotDecrementPossible(rpsCnt, lastElem, z, p, v))
                --p;
        }
        
        b = p;
        
        while ((b < lastCol) && ((z[b] == z[b + 1]) ||
               (z[b + 1] > z[b] && (rpsCnt[z[b + 1] - 1] || rpsCnt[z[b + 1]])))) {
            ++b;
        }
        
        while (BndDecrementPossible(rpsCnt, z, b))
            --b;
        
        e = b - 1;
        
        while (EdgeIncrementPossible(rpsCnt, z, e, b))
            --e;
    }
    
    PopulateVec(r, v, z, count, maxRows, isComb, partitionsVec);
    return 1;
}

int PartitionsRep(int r, const std::vector<int64_t> &v, int64_t target, int lastElem,
                  int lastCol, int maxRows, bool isComb, std::vector<int64_t> &partitionsVec) {
    
    std::vector<int> z(r);
    std::vector<int> trivVec;
    
    if (!GetFirstCombo(r, v, true, false, z, trivVec, target, trivVec, lastCol, lastElem))
        return 0;
    
    // smallest index such that z[boundary] == currMax
    int boundary = lastCol;
    
    while (boundary > 1 && z[boundary - 1] == z[lastCol])
        --boundary;
    
    // pivot is the greatest index that can be incremented.
    // We know that if z[boundary] < lastElem ==>> pivot = lastCol
    int pivot = (z[boundary] < lastElem) ? lastCol : boundary - 1;
    
    // edge is the greatest index such that z[boundary] - z[edge] >= 2
    // This is the index that will be be used as a starting point
    // to determine the next combination that meets the criteria
    int edge = boundary - 1;
    int edgeTest = z[boundary] - 2;
    
    while (edge && edgeTest < z[edge])
        --edge;
    
    int count = 0;
    
    while ((edge >= 0) && (z[boundary] - z[edge] >= 2)) {
        
        PopulateVec(r, v, z, count, maxRows, isComb, partitionsVec);
        
        if (count >= maxRows)
            break;
        
        int vertex = (z[boundary] - z[edge] == 2) ? boundary : edge + 1;
        
        ++z[edge];
        --z[vertex];
        
        if (vertex == boundary) {
            if (boundary < lastCol)
                ++boundary;
            
            const int currMax = z[boundary];
            
            while (boundary > 1 && z[boundary - 1] == currMax)
                --boundary;
            
            pivot = (z[boundary] < lastElem) ? lastCol : boundary - 1;
            
        } else if (z[vertex] == z[edge]) {
            ++vertex;
        }
        
        while (vertex < pivot) {
            const int distVert = z[vertex] - z[edge];
            const int distPivot = lastElem - z[pivot];
            
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
        
        boundary = pivot;
        
        if (boundary < lastCol && z[boundary] < z[boundary + 1])
            ++boundary;
        
        const int currMax = z[boundary];
        
        while (boundary > 1 && z[boundary - 1] == currMax)
            --boundary;
        
        edge = boundary - 1;
        edgeTest = z[boundary] - 2;
        
        while (edge && edgeTest < z[edge])
            --edge;
    }
    
    PopulateVec(r, v, z, count, maxRows, isComb, partitionsVec);
    return 1;
}

int PartitionsDistinct(int r, const std::vector<int64_t> &v, int64_t target, int lastElem,
                       int lastCol, int maxRows, bool isComb, std::vector<int64_t> &partitionsVec) {
    
    std::vector<int> z(r);
    std::vector<int> trivVec;
    
    if (!GetFirstCombo(r, v, false, false, z, trivVec, target, trivVec, lastCol, lastElem))
        return 0;
    
    int indexRows = isComb ? 0 : static_cast<int>(NumPermsNoRep(r, lastCol));
    auto indexMatrix = FromCpp14::make_unique<int[]>(indexRows * r);
    
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
    
    // Largest index such that z[boundary] - z[boundary - 1] > 1
    // This is the index that will be decremented at the same
    // time edge is incremented. See below.
    int boundary = lastCol;
    
    while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2)
        --boundary;
    
    // pivot is the greatest index that can be incremented...
    // Either z[pivot + 1] - z[pivot] > 1 or if z[lastCol] < lastElem
    // pivot = lastCol since incrementing z[lastCol] is possible
    int pivot = (z[lastCol] < lastElem) ? lastCol : boundary - 1;
    
    // edge is the greatest index such that when incremented
    // the result will be at least one less than its neighbor
    // even if its neighbor is decremented
    int edge = boundary - 1;
    int tarDiff = 3;
    
    while (edge && (z[boundary] - z[edge]) < tarDiff) {
        --edge;
        ++tarDiff;
    }
    
    int count = 0;
    
    while (edge >= 0 && (z[boundary] - z[edge]) >= tarDiff) {
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
        
        if (vertex == boundary) {
            if (boundary < lastCol)
                ++boundary;
            
            while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2)
                --boundary;
            
            pivot = (z[lastCol] < lastElem) ? lastCol : boundary - 1;
        }
        
        if (vertex < boundary) {
            if (z[vertex] - z[vertex - 1] == 1)
                ++vertex;
            
            while (vertex < pivot) {
                --z[vertex];
                ++z[pivot];
                
                if (z[vertex] - z[vertex - 1] == 1)
                    ++vertex;
                
                if (z[pivot] == lastElem || (pivot < lastCol && z[pivot + 1] - z[pivot] == 1))
                    --pivot;
            }
            
            boundary = pivot;
            
            if (boundary < lastCol && z[boundary + 1] - z[boundary] > 1)
                ++boundary;
        }
        
        edge = boundary - 1;
        tarDiff = 3;
        
        while (edge && (z[boundary] - z[edge]) < tarDiff) {
            --edge;
            ++tarDiff;
        }
    }
    
    PopulateVec(r, v, z, count, maxRows, isComb, partitionsVec);
    return 1;
}
