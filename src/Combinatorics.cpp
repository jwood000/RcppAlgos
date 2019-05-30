#include "Combinatorics.h"
#include "CleanConvert.h"
#include "RMatrix.h"
#include <RcppThread.h>

// [[Rcpp::export]]
unsigned long int cpp11GetNumThreads() {
    return std::thread::hardware_concurrency();
}

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
                           bool IsRep, typeVector myLim, double tolerance) {
    
    std::vector<int> valInd(m);
    const unsigned long int lastElem = m - 1;
    const typeVector zero = static_cast<typeVector>(0);
    
    valsData<typeVector> myVecs = valsExtreme(n, m, v, IsRep);
    std::vector<typeVector> vals = myVecs.v;
    std::vector<std::size_t> myIndex = myVecs.ind;
    const unsigned long int valSize = vals.size();
    
    // if the largest sum is smaller than myLim, there is no solution
    typeVector testMax = std::accumulate(vals.cend() - m, vals.cend(), zero);
    const std::vector<int> noSoln(1, 0);
    
    if (testMax <= myLim)  {
        if (testMax == myLim) {
            std::vector<int> maxSoln;
            
            for (auto it = myIndex.cend() - m; it < myIndex.cend(); ++it)
                maxSoln.push_back(*it);
            
            return maxSoln;
        }
        
        return noSoln;
    }
    
    // if the smallest sum is greater than myLim, there is no solution
    typeVector testMin = std::accumulate(vals.cbegin(), vals.cbegin() + m, zero);
    
    if (testMin >= myLim)  {
        if (testMin == myLim) {
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
    typeVector dist = myLim - (partial + vals[mid]);
    
    int lowBnd = (dist > 0) ? mid : 0;
    int uppBnd = (dist > 0) ? valSize - nextPos : mid;
    int ind = 0;
    std::size_t i = 0, prevInd = 0;
    
    for (; i <= lastElem; ++i) {
        while ((uppBnd - lowBnd) > 1 && dist != 0) {
            mid = (uppBnd - lowBnd) / 2;
            ind = lowBnd + mid;
            dist = myLim - (partial + vals[ind]);
            
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
            dist = myLim - (partial + vals[ind]);
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
            dist = myLim - (partial + vals[mid]);
        }
    }
    
    std::vector<int> res(m);
    std::vector<typeVector> check(m);
    
    for (std::size_t i = 0; i <= lastElem; ++i) {
        res[i] = myIndex[valInd[i]];
        check[i] = vals[valInd[i]];
    }
    
    // The algorithm above finds the first possible sum that equals
    // myLim. If there is no combination of elements from v that sum
    // to myLim, the algo returns the combination such that its sum
    // is closest to myLim and greater than myLim
    typeVector finalCheck = std::accumulate(check.cbegin(), check.cend(), zero);
    
    if (!std::is_integral<typeVector>::value) {
        if (std::abs(finalCheck - myLim) > tolerance)
            return noSoln;
    } else if (finalCheck != myLim) {
        return noSoln;
    }
    
    return res;
}

template <typename typeRcpp, typename typeVector>
typeRcpp GeneralPartitions(int n, int r, std::vector<typeVector> &v, bool isRep, typeVector lim,
                           const int numRows, bool isComb, bool xtraCol, bool bUserRows, double tol) {
    
    int count = 0;
    std::vector<typeVector> partitionsVec;
    
    std::sort(v.begin(), v.end());
    std::vector<int> z = findStart(n, r, v, isRep, lim, tol);
    
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
        partitionsVec.reserve(numRows * r);
    
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
                
                if ((numIter + count) > numRows)
                    numIter = numRows - count;
                
                for (int i = 0; i < numIter; ++i) {
                    for (int k = 0; k < r; ++k)
                        partitionsVec.push_back(v[z[k]]);
                    
                    std::next_permutation(z.begin(), z.end());
                }
                
                count += numIter;
            }
            
            if (count >= numRows)
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
        
        if (count < numRows) {
            if (isComb) {
                for (int k = 0; k < r; ++k)
                    partitionsVec.push_back(v[z[k]]);

                ++count;
            } else {
                numIter = static_cast<int>(NumPermsWithRep(z));

                if ((numIter + count) > numRows)
                    numIter = numRows - count;

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
                if (indexRows + count > numRows)
                    indexRows = numRows - count;
                
                for (int j = 0, myRow = 0; j < indexRows; ++j, myRow += r)
                    for (int k = 0; k < r; ++k)
                        partitionsVec.push_back(v[z[indexMatrix[myRow + k]]]);
                
                count += indexRows;
            }
            
            if (count >= numRows)
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
        
        if (count < numRows) {
            if (isComb) {
                for (int k = 0; k < r; ++k)
                    partitionsVec.push_back(v[z[k]]);
                
                ++count;
            } else {
                if (indexRows + count > numRows)
                    indexRows = numRows - count;
                
                for (int j = 0, myRow = 0; j < indexRows; ++j, myRow += r)
                    for (int k = 0; k < r; ++k)
                        partitionsVec.push_back(v[z[indexMatrix[myRow + k]]]);
                
                count += indexRows;
            }
        }
    }
    
    unsigned long int nCols = (xtraCol) ? r + 1 : r;
    typeRcpp partitionsMatrix = Rcpp::no_init_matrix(count, nCols);

    for (int i = 0, k = 0; i < count; ++i)
        for (int j = 0; j < r; ++j, ++k)
            partitionsMatrix(i, j) = partitionsVec[k];

    if (xtraCol)
        for (int i = 0; i < count; ++i)
            partitionsMatrix(i, r) = lim;
    
    return partitionsMatrix;
}

// This function applys a constraint function to a vector v with respect
// to a constraint value "lim". The main idea is that combinations are added
// successively, until a particular combination exceeds the given constraint
// value for a given constraint function. After this point, we can safely skip
// several combinations knowing that they will exceed the given constraint value.

template <typename typeRcpp, typename typeVector>
typeRcpp CombinatoricsConstraints(int n, int r, std::vector<typeVector> &v, bool isRep, std::string myFun,
                                  std::vector<std::string> comparison, std::vector<typeVector> lim, int numRows,
                                  bool isComb, bool xtraCol, std::vector<int> &Reps, bool isMult, double tol) {
    
    // myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max";
    // The comparison vector is a comparison operator: 
    //                             "<", "<=", ">", ">=", "==", ">,<", ">=,<", ">,<=", ">=,<=";
    
    typeVector testVal;
    int count = 0;
    const int numCols = xtraCol ? (r + 1) : r;
    unsigned long int uR = r;
    typeRcpp combinatoricsMatrix = Rcpp::no_init_matrix(numRows, numCols);
    
    Rcpp::XPtr<funcPtr<typeVector>> xpFun = putFunPtrInXPtr<typeVector>(myFun);
    funcPtr<typeVector> constraintFun = *xpFun;
    
    // We first check if we are getting double precision.
    // If so, for the non-strict inequalities, we have
    // to alter the limit by epsilon:
    //
    //           x <= y   --->>>   x <= y + e
    //           x >= y   --->>>   x >= y - e
    //
    // Equality is a bit tricky as we need to check
    // whether the absolute value of the difference is
    // less than epsilon. As a result, we can't alter
    // the limit with one alteration. Observe:
    //
    //   x == y  --->>>  |x - y| <= e , which gives:
    //
    //             - e <= x - y <= e
    //
    //         1.     x >= y - e
    //         2.     x <= y + e
    //
    // As a result, we must define a specialized equality
    // check for double precision. It is 'equalDbl' and
    // can be found in ConstraintsUtils.h
    
    if (!std::is_integral<typeVector>::value) {
        if (comparison[0] == "==") {
            lim.push_back(lim[0] - tol);
            lim[0] += tol;
        }
        
        if (comparison[0] == "<=") {
            lim[0] += tol;
        } else if (comparison[0] == ">=") {
            lim[0] -= tol;
        }
        
        if (comparison.size() > 1) {
            if (comparison[1] == "<=") {
                lim[1] += tol;
            } else if (comparison[1] == ">=") {
                lim[1] -= tol;
            }
        }
    }
    
    for (std::size_t nC = 0; nC < comparison.size(); ++nC) {
        
        Rcpp::XPtr<compPtr<typeVector>> xpCompOne = putCompPtrInXPtr<typeVector>(comparison[nC]);
        compPtr<typeVector> comparisonFunOne = *xpCompOne;
        
        Rcpp::XPtr<compPtr<typeVector>> xpCompTwo = xpCompOne;
        compPtr<typeVector> comparisonFunTwo;
    
        if (comparison[nC] == ">" || comparison[nC] == ">=") {
            if (isMult) {
                for (int i = 0; i < (n - 1); ++i) {
                    for (int j = (i + 1); j < n; ++j) {
                        if (v[i] < v[j]) {
                            std::swap(v[i], v[j]);
                            std::swap(Reps[i], Reps[j]);
                        }
                    }
                }
            } else {
                std::sort(v.begin(), v.end(), std::greater<double>());
            }
            comparisonFunTwo = *xpCompOne;
        } else {
            if (isMult) {
                for (int i = 0; i < (n-1); ++i) {
                    for (int j = (i+1); j < n; ++j) {
                        if (v[i] > v[j]) {
                            std::swap(v[i], v[j]);
                            std::swap(Reps[i], Reps[j]);
                        }
                    }
                }
            } else {
                std::sort(v.begin(), v.end());
            }
            
            std::vector<std::string>::const_iterator itComp = std::find(compSpecial.cbegin(), 
                                                                        compSpecial.cend(), 
                                                                        comparison[nC]);
            if (itComp != compSpecial.end()) {
                int myIndex = std::distance(compSpecial.cbegin(), itComp);
                Rcpp::XPtr<compPtr<typeVector>> xpCompThree = putCompPtrInXPtr<typeVector>(compHelper[myIndex]);
                comparisonFunTwo = *xpCompThree;
            } else {
                comparisonFunTwo = *xpCompOne;
            }
        }
        
        std::vector<int> z, zCheck;
        std::vector<typeVector> testVec(r);
        bool t_1, t_2, t = true, keepGoing = true;
        int numIter, myStart, maxZ = n - 1;
        const int r1 = r - 1;
        const int r2 = r - 2;
        
        if (isMult) {
            int zExpSize = std::accumulate(Reps.cbegin(), Reps.cend(), 0);
            std::vector<int> zExpand, zIndex, zGroup(r), zPerm(r);
            
            for (int i = 0, k = 0; i < n; ++i) {
                zIndex.push_back(k);
                
                for (int j = 0; j < Reps[i]; ++j, ++k)
                    zExpand.push_back(i);
            }
            
            for (int i = 0; i < r; ++i)
                z.push_back(zExpand[i]);
            
            while (keepGoing) {
                
                t_2 = true;
                for (int i = 0; i < r; ++i)
                    testVec[i] = v[zExpand[zIndex[z[i]]]];
                
                testVal = constraintFun(testVec, uR);
                t = comparisonFunTwo(testVal, lim);
                
                while (t && t_2 && keepGoing) {
                    
                    testVal = constraintFun(testVec, uR);
                    t_1 = comparisonFunOne(testVal, lim);
                    
                    if (t_1) {
                        myStart = count;
                        
                        if (isComb) {
                            for (int k = 0; k < r; ++k)
                                combinatoricsMatrix(count, k) = v[zExpand[zIndex[z[k]]]];
                            
                            ++count;
                        } else {
                            for (int k = 0; k < r; ++k)
                                zPerm[k] = zExpand[zIndex[z[k]]];
                            
                            numIter = static_cast<int>(NumPermsWithRep(zPerm));
                            
                            if ((numIter + count) > numRows)
                                numIter = numRows - count;
                            
                            for (int i = 0; i < numIter; ++i) {
                                for (int k = 0; k < r; ++k)
                                    combinatoricsMatrix(count, k) = v[zPerm[k]];
                                
                                std::next_permutation(zPerm.begin(), zPerm.end());
                            }
                            
                            count += numIter;
                        }
                        
                        if (xtraCol)
                            for (int i = myStart; i < count; ++i)
                                combinatoricsMatrix(i, r) = testVal;
                    }
                    
                    keepGoing = (count < numRows);
                    t_2 = (z[r1] != maxZ);
                    
                    if (t_2) {
                        ++z[r1];
                        testVec[r1] = v[zExpand[zIndex[z[r1]]]];
                        testVal = constraintFun(testVec, uR);
                        t = comparisonFunTwo(testVal, lim);
                    }
                }
                
                if (keepGoing) {
                    zCheck = z;
                    for (int i = r2; i >= 0; --i) {
                        if (zExpand[zIndex[z[i]]] != zExpand[zExpSize - r + i]) {
                            ++z[i];
                            testVec[i] = v[zExpand[zIndex[z[i]]]];
                            zGroup[i] = zIndex[z[i]];
                            
                            for (int k = (i+1); k < r; ++k) {
                                zGroup[k] = zGroup[k-1] + 1;
                                z[k] = zExpand[zGroup[k]];
                                testVec[k] = v[zExpand[zIndex[z[k]]]];
                            }
                            
                            testVal = constraintFun(testVec, uR);
                            t = comparisonFunTwo(testVal, lim);
                            if (t) {break;}
                        }
                    }
                    
                    if (!t || zCheck == z) {keepGoing = false;}
                }
            }
            
        } else if (isRep) {
            
            v.erase(std::unique(v.begin(), v.end()), v.end());
            z.assign(r, 0);
            maxZ = static_cast<int>(v.size()) - 1;
            
            while (keepGoing) {
                
                t_2 = true;
                for (int i = 0; i < r; ++i)
                    testVec[i] = v[z[i]];
                
                testVal = constraintFun(testVec, uR);
                t = comparisonFunTwo(testVal, lim);
                
                while (t && t_2 && keepGoing) {
                    
                    testVal = constraintFun(testVec, uR);
                    t_1 = comparisonFunOne(testVal, lim);
                    
                    if (t_1) {
                        myStart = count;
                        
                        if (isComb) {
                            for (int k = 0; k < r; ++k)
                                combinatoricsMatrix(count, k) = v[z[k]];
                            
                            ++count;
                        } else {
                            numIter = static_cast<int>(NumPermsWithRep(z));
                            
                            if ((numIter + count) > numRows)
                                numIter = numRows - count;
                            
                            for (int i = 0; i < numIter; ++i) {
                                for (int k = 0; k < r; ++k)
                                    combinatoricsMatrix(count, k) = v[z[k]];
                                
                                std::next_permutation(z.begin(), z.end());
                            }
                            
                            count += numIter;
                        }
                        
                        if (xtraCol)
                            for (int i = myStart; i < count; ++i)
                                combinatoricsMatrix(i, r) = testVal;
                        
                        keepGoing = (count < numRows);
                    }
                    
                    t_2 = (z[r1] != maxZ);
                    
                    if (t_2) {
                        ++z[r1];
                        testVec[r1] = v[z[r1]];
                        testVal = constraintFun(testVec, uR);
                        t = comparisonFunTwo(testVal, lim);
                    }
                }
                
                if (keepGoing) {
                    zCheck = z;
                    for (int i = r2; i >= 0; --i) {
                        if (z[i] != maxZ) {
                            ++z[i];
                            testVec[i] = v[z[i]];
                            
                            for (int k = (i+1); k < r; ++k) {
                                z[k] = z[k-1];
                                testVec[k] = v[z[k]];
                            }
                            
                            testVal = constraintFun(testVec, uR);
                            t = comparisonFunTwo(testVal, lim);
                            if (t) {break;}
                        }
                    }
                    
                    if (!t || zCheck == z) {keepGoing = false;}
                }
            }
            
        } else {
            
            for (int i = 0; i < r; ++i)
                z.push_back(i);
            
            const int nMinusR = (n - r);
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
    
            while (keepGoing) {
    
                t_2 = true;
                for (int i = 0; i < r; ++i)
                    testVec[i] = v[z[i]];
    
                testVal = constraintFun(testVec, uR);
                t = comparisonFunTwo(testVal, lim);
    
                while (t && t_2 && keepGoing) {
    
                    testVal = constraintFun(testVec, uR);
                    t_1 = comparisonFunOne(testVal, lim);
    
                    if (t_1) {
                        myStart = count;
                        
                        if (isComb) {
                            for (int k = 0; k < r; ++k)
                                combinatoricsMatrix(count, k) = v[z[k]];
    
                            ++count;
                        } else {
                            if (indexRows + count > numRows)
                                indexRows = numRows - count;
                            
                            for (int j = 0, myRow = 0; j < indexRows; ++j, myRow += r)
                                for (int k = 0; k < r; ++k)
                                    combinatoricsMatrix(count, k) = v[z[indexMatrix[myRow + k]]];
                            
                            count += indexRows;
                        }
                        
                        if (xtraCol)
                            for (int i = myStart; i < count; ++i)
                                combinatoricsMatrix(i, r) = testVal;
                        
                        keepGoing = (count < numRows);
                    }
    
                    t_2 = (z[r1] != maxZ);
    
                    if (t_2) {
                        ++z[r1];
                        testVec[r1] = v[z[r1]];
                        testVal = constraintFun(testVec, uR);
                        t = comparisonFunTwo(testVal, lim);
                    }
                }
    
                if (keepGoing) {
                    zCheck = z;
                    for (int i = r2; i >= 0; --i) {
                        if (z[i] != (nMinusR + i)) {
                            ++z[i];
                            testVec[i] = v[z[i]];
                            
                            for (int k = (i+1); k < r; ++k) {
                                z[k] = z[k - 1] + 1;
                                testVec[k] = v[z[k]];
                            }
                            
                            testVal = constraintFun(testVec, uR);
                            t = comparisonFunTwo(testVal, lim);
                            if (t) {break;}
                        }
                    }
                    
                    if (!t || zCheck == z) {keepGoing = false;}
                }
            }
        }
        
        lim.erase(lim.begin());
    }
       
    return SubMat(combinatoricsMatrix, count);
}

// [[Rcpp::export]]
SEXP CombinatoricsRcpp(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow,
                       SEXP Rhigh, SEXP f1, SEXP f2, SEXP Rlim, bool IsComb, 
                       SEXP RKeepRes, bool IsFactor, bool IsCount, SEXP stdFun, 
                       SEXP myEnv, SEXP Rparallel, SEXP RNumThreads, int maxThreads,
                       SEXP Rtolerance) {
    
    int n, m1, m2, m = 0, lenFreqs = 0, nRows = 0;
    bool IsLogical, IsCharacter, IsMultiset, IsInteger;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqsExpanded;
    Rcpp::CharacterVector rcppChar;
    
    bool keepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRepetition = CleanConvert::convertLogical(RisRep, "repetition");
    
    switch(TYPEOF(Rv)) {
        case LGLSXP: {
            IsLogical = true;
            keepRes = IsInteger = IsCharacter = false;
            break;
        }
        case INTSXP: {
            IsInteger = true;
            IsLogical = IsCharacter = false;
            break;
        }
        case REALSXP: {
            IsLogical = IsInteger = IsCharacter = false;
            break;
        }
        case STRSXP: {
            IsCharacter = true;
            Parallel = keepRes = IsLogical = IsInteger = false;
            break;
        }
        default: {
            Rcpp::stop("Only integers, numerical, character, and factor classes are supported for v");   
        }
    }
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        myReps.push_back(1);
    } else {
        IsMultiset = true;
        IsRepetition = false;
        CleanConvert::convertVector(RFreqs, myReps, "freqs");
        lenFreqs = static_cast<int>(myReps.size());
        
        for (int i = 0; i < lenFreqs; ++i)
            for (int j = 0; j < myReps[i]; ++j)
                freqsExpanded.push_back(i);
    }
    
    if (Rf_isNull(Rm)) {
        if (IsMultiset) {
            m = freqsExpanded.size();
        } else {
            Rcpp::stop("m and freqs cannot both be NULL");
        }
    } else {
        if (Rf_length(Rm) > 1)
            Rcpp::stop("length of m must be 1");
        
        CleanConvert::convertPrimitive(Rm, m, "m");
    }
    
    std::vector<double> rowVec(m);
    
    if (IsCharacter) {
        rcppChar = Rcpp::as<Rcpp::CharacterVector>(Rv);
        n = rcppChar.size();
    } else if (IsLogical) {
        vInt = Rcpp::as<std::vector<int>>(Rv);
        n = vInt.size();
    } else {
        if (Rf_length(Rv) == 1) {
            int seqEnd;             // numOnly = true, checkWhole = true, negPoss = true
            CleanConvert::convertPrimitive(Rv, seqEnd, "If v is not a character and of length 1, it", true, true, true);
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            Rcpp::IntegerVector vTemp = Rcpp::seq(m1, m2);
            IsInteger = true;
            vNum = Rcpp::as<std::vector<double>>(vTemp);
        } else {
            vNum = Rcpp::as<std::vector<double>>(Rv);
        }
        
        n = vNum.size();
    }
    
    if (IsInteger) {
        for (int i = 0; i < n && IsInteger; ++i)
            if (Rcpp::NumericVector::is_na(vNum[i]))
                IsInteger = false;
        
        if (IsInteger)
            vInt.assign(vNum.cbegin(), vNum.cend());
    }
        
    bool IsConstrained;
    if (IsFactor)
        keepRes = IsConstrained = IsCharacter = IsInteger = false;
    
    if (IsLogical || IsCharacter || Rf_isNull(f1) || Rf_isNull(f2) || Rf_isNull(Rlim)) {
        IsConstrained = false;
    } else {
        if (!Rf_isString(f1))
            Rcpp::stop("constraintFun must be passed as a character");
        
        if (!Rf_isString(f2))
            Rcpp::stop("comparisonFun must be passed as a character");
        
        IsConstrained = true;
    }
    
    if (IsConstrained) {
        for (int i = (vNum.size() - 1); i >= 0; --i)
            if (Rcpp::NumericVector::is_na(vNum[i]))
                vNum.erase(vNum.begin() + i);
            
        n = vNum.size();
    }

    double computedRows = 0;
    
    if (IsMultiset) {
        if (n != lenFreqs)
            Rcpp::stop("the length of freqs must equal the length of v");
        
        if (m > static_cast<int>(freqsExpanded.size()))
            m = freqsExpanded.size();
        
        if (IsComb) {
            computedRows = MultisetCombRowNum(n, m, myReps);
        } else {
            if (Rf_isNull(Rm) || m == static_cast<int>(freqsExpanded.size()))
                computedRows = NumPermsWithRep(freqsExpanded);
            else
                computedRows = MultisetPermRowNum(n, m, myReps);
        }
    } else {
        if (IsRepetition) {
            if (IsComb)
                computedRows = NumCombsWithRep(n, m);
            else
                computedRows = std::pow((double) n, (double) m);
        } else {
            if (m > n)
                Rcpp::stop("m must be less than or equal to the length of v");
            
            if (IsComb)
                computedRows = nChooseK(n, m);
            else
                computedRows = NumPermsNoRep(n, m);
        }
    }
    
    bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);
    
    if (IsGmp) {
        if (IsMultiset) {
            if (IsComb) {
                MultisetCombRowNumGmp(computedRowMpz, n, m, myReps);
            } else {
                if (Rf_isNull(Rm) || m == static_cast<int>(freqsExpanded.size()))
                    NumPermsWithRepGmp(computedRowMpz, freqsExpanded);
                else
                    MultisetPermRowNumGmp(computedRowMpz, n, m, myReps);
            }
        } else {
            if (IsRepetition) {
                if (IsComb)
                    NumCombsWithRepGmp(computedRowMpz, n, m);
                else
                    mpz_ui_pow_ui(computedRowMpz, n, m);
            } else {
                if (IsComb)
                    nChooseKGmp(computedRowMpz, n, m);
                else
                    NumPermsNoRepGmp(computedRowMpz, n, m);
            }
        }
    }

    double lower = 0, upper = 0;
    bool bLower = false, bUpper = false;
    mpz_t lowerMpz[1], upperMpz[1];
    mpz_init(lowerMpz[0]); mpz_init(upperMpz[0]);
    mpz_set_ui(lowerMpz[0], 0); mpz_set_ui(upperMpz[0], 0);
    
    if (!IsCount) {
        if (IsGmp) {
            if (!Rf_isNull(Rlow)) {
                bLower = true;
                createMPZArray(Rlow, lowerMpz, 1, "lower");
                mpz_sub_ui(lowerMpz[0], lowerMpz[0], 1);
            }
            
            if (!Rf_isNull(Rhigh)) {
                bUpper = true;
                createMPZArray(Rhigh, upperMpz, 1, "upper");
            }
        } else { 
            if (!Rf_isNull(Rlow)) {
                bLower = true;                        // numOnly = false
                CleanConvert::convertPrimitive(Rlow, lower, "lower", false);
                --lower;
            }
            
            if (!Rf_isNull(Rhigh)) {
                bUpper = true;                           // numOnly = false
                CleanConvert::convertPrimitive(Rhigh, upper, "upper", false);
            }
        }
    }
    
    if (IsGmp) {
        if (mpz_cmp(lowerMpz[0], computedRowMpz) >= 0 || mpz_cmp(upperMpz[0], computedRowMpz) > 0)
            Rcpp::stop("bounds cannot exceed the maximum number of possible results");
    } else {
        if (lower >= computedRows || upper > computedRows)
            Rcpp::stop("bounds cannot exceed the maximum number of possible results");
    }
    
    if (IsCount) {
        if (IsGmp) {
            unsigned long int sizeNum, size = sizeof(int);
            unsigned long int numb = 8 * sizeof(int);
            sizeNum = sizeof(int) * (2 + (mpz_sizeinbase(computedRowMpz, 2) + numb - 1) / numb);
            size += sizeNum;
            
            SEXP ansPos = PROTECT(Rf_allocVector(RAWSXP, size));
            char* rPos = (char*)(RAW(ansPos));
            ((int*)(rPos))[0] = 1; // first int is vector-size-header
            
            // current position in rPos[] (starting after vector-size-header)
            unsigned long int posPos = sizeof(int);
            posPos += myRaw(&rPos[posPos], computedRowMpz, sizeNum);
            
            Rf_setAttrib(ansPos, R_ClassSymbol, Rf_mkString("bigz"));
            UNPROTECT(1);
            return(ansPos);
        } else {
            if (computedRows > std::numeric_limits<int>::max())
                return Rcpp::wrap(computedRows);
            else
                return Rcpp::wrap(static_cast<int>(computedRows));
        }
    }
    
    std::vector<int> startZ(m);
    bool permNonTriv = false;
    double dblLower = lower;
    if (!IsGmp) mpz_set_d(lowerMpz[0], dblLower);
    
    if (bLower && mpz_cmp_ui(lowerMpz[0], 0) > 0) {
        if (IsComb) {
            if (IsGmp)
                startZ = nthCombinationGmp(n, m, lowerMpz[0], IsRepetition, IsMultiset, myReps);
            else
                startZ = nthCombination(n, m, lower, IsRepetition, IsMultiset, myReps);
        } else {
            permNonTriv = true;
            if (IsGmp)
                startZ = nthPermutationGmp(n, m, lowerMpz[0], IsRepetition, IsMultiset, myReps, freqsExpanded, true);
            else
                startZ = nthPermutation(n, m, lower, IsRepetition, IsMultiset, myReps, freqsExpanded, true);
        }
    } else {
        if (IsComb) {
            if (IsMultiset)
                startZ.assign(freqsExpanded.cbegin(), freqsExpanded.cbegin() + m);
            else if (IsRepetition)
                std::fill(startZ.begin(), startZ.end(), 0);
            else
                std::iota(startZ.begin(), startZ.end(), 0);
        } else {
            if (IsMultiset) {
                startZ = freqsExpanded;
            } else if (IsRepetition) {
                std::fill(startZ.begin(), startZ.end(), 0);
            } else {
                startZ.resize(n);
                std::iota(startZ.begin(), startZ.end(), 0);
            }
        }
    }
    
    double userNumRows = 0;
    
    if (IsGmp) {
        mpz_t testBound; mpz_init(testBound);
        if (bLower && bUpper) {
            mpz_sub(testBound, upperMpz[0], lowerMpz[0]);
            mpz_abs(testBound, testBound);
            if (mpz_cmp_ui(testBound, std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(testBound);
        } else if (bUpper) {
            permNonTriv = true;
            if (mpz_cmp_d(upperMpz[0], std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
                
            userNumRows = mpz_get_d(upperMpz[0]);
        } else if (bLower) {
            mpz_sub(testBound, computedRowMpz, lowerMpz[0]);
            mpz_abs(testBound, testBound);
            if (mpz_cmp_d(testBound, std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(testBound);
        }
        mpz_clear(testBound);
    } else {
        if (bLower && bUpper)
            userNumRows = upper - lower;
        else if (bUpper)
            userNumRows = upper;
        else if (bLower)
            userNumRows = computedRows - lower;
    }
    
    if (userNumRows == 0) {
        if (bLower && bUpper) {
            // Since lower is decremented and upper isn't upper - lower = 0 means
            // that lower is one larger than upper as put in by the user
            
            Rcpp::stop("The number of rows must be positive. Either the lowerBound "
                           "exceeds the maximum number of possible results or the "
                           "lowerBound is greater than the upperBound.");
        } else {
            if (computedRows > std::numeric_limits<int>::max())
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = computedRows;
            nRows = static_cast<int>(computedRows);
        }
    } else if (userNumRows < 0) {
        Rcpp::stop("The number of rows must be positive. Either the lowerBound "
                  "exceeds the maximum number of possible results or the "
                  "lowerBound is greater than the upperBound.");
    } else if (userNumRows > std::numeric_limits<int>::max()) {
        Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
    } else {
        nRows = static_cast<int>(userNumRows);
    }
    
    unsigned long int uM = m;
    int nThreads = 1;
    
    // Determined empirically. Setting up threads can be expensive,
    // so we set the cutoff below to ensure threads aren't spawned
    // unnecessarily. We also protect users with fewer than 2 threads
    if ((nRows < 20000) || (maxThreads < 2)) {
        Parallel = false;
    } else if (!Rf_isNull(RNumThreads)) {
        int userThreads = 1;
        if (!Rf_isNull(RNumThreads))
            CleanConvert::convertPrimitive(RNumThreads, userThreads, "nThreads");
        
        if (userThreads > maxThreads)
            userThreads = maxThreads;
        
        // Ensure that each thread has at least 10000
        if ((nRows / userThreads) < 10000)
            userThreads = nRows / 10000;

        if (userThreads > 1 && !IsCharacter) {
            Parallel = true;
            nThreads = userThreads;
        }
    } else if (Parallel) {
        nThreads = (maxThreads > 2) ? (maxThreads - 1) : 2;
        
        // Ensure that each thread has at least 10000
        if ((nRows / nThreads) < 10000)
            nThreads = nRows / 10000;
    }
    
    if (IsConstrained) {
        std::vector<double> myLim;       // numOnly = true, checkWhole = false, negPoss = true
        CleanConvert::convertVector(Rlim, myLim, "limitConstraints", true, false, true);
        
        if (myLim.size() > 2)
            Rcpp::stop("there cannot be more than 2 limitConstraints");
        else if (myLim.size() == 2 && myLim[0] == myLim[1])
            Rcpp::stop("The limitConstraints must be different");
        
        std::string mainFun = Rcpp::as<std::string>(f1);
        if (mainFun != "prod" && mainFun != "sum" && mainFun != "mean"
                && mainFun != "max" && mainFun != "min") {
            Rcpp::stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
        }

        std::vector<std::string>::const_iterator itComp;
        std::vector<std::string> compFunVec = Rcpp::as<std::vector<std::string>>(f2);
        
        if (compFunVec.size() > 2)
            Rcpp::stop("there cannot be more than 2 comparison operators");
        
        for (std::size_t i = 0; i < compFunVec.size(); ++i) {
            itComp = std::find(compForms.cbegin(), compForms.cend(), compFunVec[i]);
            
            if (itComp == compForms.end())
                Rcpp::stop("comparison operators must be one of the following: '>', '>=', '<', '<=', or '=='");
                
            int myIndex = std::distance(compForms.cbegin(), itComp);
            
            // The first 5 are "standard" whereas the 6th and 7th
            // are written with the equality first. Converting
            // them here makes it easier to deal with later.
            if (myIndex > 4)
                myIndex -= 3;
            
            compFunVec[i] = compForms[myIndex];
        }
        
        if (compFunVec.size() == 2) {
            if (myLim.size() == 1) {
                compFunVec.pop_back();
            } else {
                if (compFunVec[0] == "==" || compFunVec[1] == "==")
                    Rcpp::stop("If comparing against two limitConstraints, the "
                                "equality comparisonFun (i.e. '==') cannot be used. "
                                "Instead, use '>=' or '<='.");
                
                if (compFunVec[0].substr(0, 1) == compFunVec[1].substr(0, 1))
                    Rcpp::stop("Cannot have two 'less than' comparisonFuns or two 'greater than' "
                          "comparisonFuns (E.g. c('<', '<=') is not allowed).");
                
                bool sortNeeded = false;
                
                if (compFunVec[0].substr(0, 1) == ">" && std::min(myLim[0], myLim[1]) == myLim[0]) {
                    compFunVec[0] = compFunVec[0] + "," + compFunVec[1];
                    sortNeeded = true;
                } else if (compFunVec[0].substr(0, 1) == "<" && std::max(myLim[0], myLim[1]) == myLim[0]) {
                    compFunVec[0] = compFunVec[1] + "," + compFunVec[0];
                    sortNeeded = true;
                }
                
                if (sortNeeded) {
                    if (std::max(myLim[0], myLim[1]) == myLim[1]) {
                        double temp = myLim[0];
                        myLim[0] = myLim[1];
                        myLim[1] = temp;
                    }
                    
                    compFunVec.pop_back();
                }
            }
        } else {
            if (myLim.size() == 2)
                myLim.pop_back();
        }
        
        bool SpecialCase = false;
        
        // If bLower, the user is looking to test a particular range. Otherwise, the constraint algo
        // will simply return (upper - lower) # of combinations/permutations that meet the criteria
        if (bLower) {
            SpecialCase = true;
        } else if (mainFun == "prod") {
            for (int i = 0; i < n; ++i) {
                if (vNum[i] < 0) {
                    SpecialCase = true;
                    break;
                }
            }
        }
        
        Rcpp::XPtr<funcPtr<double>> xpFunDbl = putFunPtrInXPtr<double>(mainFun);
        funcPtr<double> myFunDbl = *xpFunDbl;
        
        IsInteger = (IsInteger) && checkIsInteger(mainFun, uM, n, rowVec, vNum, myLim, myFunDbl, true);
        // Must be defined inside IsInteger check as myLim could be outside integer data type range
        std::vector<int> limInt;
        double tolerance = 0;
        
        if (IsInteger) {
            limInt.assign(myLim.cbegin(), myLim.cend());
        } else {
            CleanConvert::convertPrimitive(Rtolerance, tolerance, "tolerance", true, false, false, true);
        }
        
        if (SpecialCase) {
            if (IsInteger) {
                return SpecCaseRet<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, nRows, keepRes, startZ, lower,
                                                        mainFun, IsMultiset, computedRows, compFunVec, limInt, IsComb,
                                                        myReps, freqsExpanded, bLower, permNonTriv, userNumRows, tolerance);
            } else {
                return SpecCaseRet<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes, startZ, lower,
                                                        mainFun, IsMultiset, computedRows, compFunVec, myLim, IsComb,
                                                        myReps, freqsExpanded, bLower, permNonTriv, userNumRows, tolerance);
            }
        }
        
        if (mainFun == "sum" && compFunVec[0] == "==" && !IsMultiset && nRows > 100 && n > 1) {
            bool PartitionCase = true;
            std::vector<double> pTest(vNum.cbegin(), vNum.cend());
            std::sort(pTest.begin(), pTest.end());
            const double tarDiff = pTest[1] - pTest[0];
            constexpr double limitTest = 2 * std::numeric_limits<double>::epsilon();
            
            if (static_cast<int64_t>(pTest[0]) == pTest[0]) {
                for (int i = 1; i < n; ++i) {
                    const double testDiff = pTest[i] - pTest[i - 1];
                    
                    if (std::abs(testDiff - tarDiff) > limitTest 
                            || static_cast<int64_t>(pTest[i]) != pTest[i]) {
                        PartitionCase = false;
                        break;
                    }
                }
                
                if (PartitionCase) {
                    bool bUserRows = bLower || bUpper;
                    
                    if (IsInteger) {
                        return GeneralPartitions<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, limInt[0],
                                                                      nRows, IsComb, keepRes, bUserRows, tolerance);
                    } else {
                        return GeneralPartitions<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, myLim[0],
                                                                      nRows, IsComb, keepRes, bUserRows, tolerance);
                    }
                }
            }
        }

        if (IsInteger) {
            return CombinatoricsConstraints<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, mainFun, compFunVec,
                                                                 limInt, nRows, IsComb, keepRes, myReps, IsMultiset, tolerance);
        }
        
        return CombinatoricsConstraints<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, mainFun, compFunVec,
                                                             myLim, nRows, IsComb, keepRes, myReps, IsMultiset, tolerance);
    } else {
        bool applyFun = !Rf_isNull(stdFun) && !IsFactor;

        if (applyFun) {
            if (!Rf_isFunction(stdFun))
                Rcpp::stop("FUN must be a function!");
            
            SEXP ans = PROTECT(Rf_allocVector(VECSXP, nRows));
            SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));
            
            if (IsCharacter) {
                ApplyFunction(n, m, rcppChar, IsRepetition, nRows, IsComb, myReps, 
                              ans, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            } else if (IsLogical || IsInteger) {
                Rcpp::IntegerVector rcppVInt(vInt.cbegin(), vInt.cend());
                ApplyFunction(n, m, rcppVInt, IsRepetition, nRows, IsComb, myReps, 
                              ans, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            } else {
                Rcpp::NumericVector rcppVNum(vNum.cbegin(), vNum.cend());
                ApplyFunction(n, m, rcppVNum, IsRepetition, nRows, IsComb, myReps, 
                              ans, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            }
            
            UNPROTECT(2);
            return ans;
        }
        
        // It is assumed that if user has constraintFun with no comparison
        // or limitConstraints and they have not explicitly set keepRes to
        // FALSE, then they simply want the constraintFun applied
        if (Rf_isNull(RKeepRes)) {
            if (Rf_isNull(f2) && Rf_isNull(Rlim) && !Rf_isNull(f1))
                keepRes = !IsLogical && !IsCharacter && !IsFactor;
        } else {
            keepRes = keepRes && !Rf_isNull(f1);
        }
        
        std::string mainFun;
        funcPtr<double> myFunDbl;
        funcPtr<int> myFunInt;
        int nCol = m;
        
        if (keepRes) {
            mainFun = Rcpp::as<std::string>(f1);
            if (mainFun != "prod" && mainFun != "sum" && mainFun != "mean"
                    && mainFun != "max" && mainFun != "min") {
                Rcpp::stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
            }
            
            Rcpp::XPtr<funcPtr<double>> xpFunDbl = putFunPtrInXPtr<double>(mainFun);
            Rcpp::XPtr<funcPtr<int>> xpFunInt = putFunPtrInXPtr<int>(mainFun);
            myFunDbl = *xpFunDbl;
            myFunInt = *xpFunInt;
            
            IsInteger = (IsInteger) && checkIsInteger(mainFun, uM, n, rowVec, vNum, vNum, myFunDbl);
            ++nCol;
        }
        
        if (Parallel) {
            permNonTriv = true;
            RcppThread::ThreadPool pool(nThreads);
            int step = 0, stepSize = nRows / nThreads;
            int nextStep = stepSize;
            
            if (IsLogical) {
                Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(nRows, nCol);
                RcppParallel::RMatrix<int> parBool(matBool);

                for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                    pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<int>, int>), 
                              n, m, vInt, IsRepetition, nextStep, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunInt, keepRes, std::ref(parBool), step);
                    
                    getStartZ(n, m, lower, stepSize, lowerMpz[0], IsRepetition,
                              IsComb, IsMultiset, IsGmp, myReps, freqsExpanded, startZ);
                }

                pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<int>, int>), 
                          n, m, vInt, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                          startZ, permNonTriv, IsMultiset, myFunInt, keepRes, std::ref(parBool), step);

                pool.join();
                return matBool;
                
            } else if (IsFactor || IsInteger) {
                Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(nRows, nCol);
                RcppParallel::RMatrix<int> parInt(matInt);
                
                for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                    pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<int>, int>), 
                              n, m, vInt, IsRepetition, nextStep, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunInt, keepRes, std::ref(parInt), step);

                    getStartZ(n, m, lower, stepSize, lowerMpz[0], IsRepetition, 
                               IsComb, IsMultiset, IsGmp, myReps, freqsExpanded, startZ);
                }
                
                pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<int>, int>),
                          n, m, vInt, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                          startZ, permNonTriv, IsMultiset, myFunInt, keepRes, std::ref(parInt), step);
                
                pool.join();
                
                if (IsFactor) {
                    Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
                    Rcpp::CharacterVector myClass = testFactor.attr("class");
                    Rcpp::CharacterVector myLevels = testFactor.attr("levels");
                    matInt.attr("class") = myClass;
                    matInt.attr("levels") = myLevels;
                }
                
                return matInt;
                
            } else {
                Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(nRows, nCol);
                RcppParallel::RMatrix<double> parNum(matNum);
                
                for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                    pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<double>, double>), 
                              n, m, vNum, IsRepetition, nextStep, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunDbl, keepRes, std::ref(parNum), step);
                    
                    getStartZ(n, m, lower, stepSize, lowerMpz[0], IsRepetition, 
                              IsComb, IsMultiset, IsGmp, myReps, freqsExpanded, startZ);
                }

                pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<double>, double>),
                          n, m, vNum, IsRepetition, nRows, IsComb, myReps, freqsExpanded, startZ,
                          permNonTriv, IsMultiset, myFunDbl, keepRes, std::ref(parNum), step);
                
                pool.join();
                return matNum;
            }
        } else {
            if (IsCharacter) {
                Rcpp::CharacterMatrix matChar = Rcpp::no_init_matrix(nRows, nCol);
                CharacterReturn(n, m, rcppChar, IsRepetition, nRows, IsComb, myReps,
                                freqsExpanded, startZ, permNonTriv, IsMultiset, keepRes, matChar, 0);
                return matChar;
                
            } else if (IsLogical) {
                Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(nRows, nCol);
                GeneralReturn(n, m, vInt, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunInt, keepRes, matBool, 0);
                return matBool;
                
            } else if (IsFactor || IsInteger) {
                Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(nRows, nCol);
                GeneralReturn(n, m, vInt, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunInt, keepRes, matInt, 0);
                
                if (IsFactor) {
                    Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
                    Rcpp::CharacterVector myClass = testFactor.attr("class");
                    Rcpp::CharacterVector myLevels = testFactor.attr("levels");
                    matInt.attr("class") = myClass;
                    matInt.attr("levels") = myLevels;
                }
                
                return matInt;
                
            } else {
                Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(nRows, nCol);
                GeneralReturn(n, m, vNum, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunDbl, keepRes, matNum, 0);
                return matNum;
            }
        }
    }
}
