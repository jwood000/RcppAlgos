#ifndef CONSTRAINTS_MASTER_H
#define CONSTRAINTS_MASTER_H

#include "Combinations.h"
#include "Permutations.h"
#include "ComboResults.h"
#include "PermuteResults.h"
#include "GeneralPartitions.h"

const std::vector<std::string> compForms = {"<", ">", "<=", ">=", "==", "=<", "=>"};
const std::vector<std::string> compSpecial = {"==", ">,<", ">=,<", ">,<=", ">=,<="};
const std::vector<std::string> compHelper = {"<=", "<", "<", "<=", "<="};

// N.B. Passing signed int to function expecting std::size_t is well defined
template <typename typeRcpp, typename typeVector>
void GeneralReturn(int n, int r, std::vector<typeVector> v, bool IsRep, int nRows, bool IsComb,
                   std::vector<int> myReps, std::vector<int> freqs, std::vector<int> z,
                   bool generalReturn, bool IsMultiset, funcPtr<typeVector> myFun,
                   bool keepRes, typeRcpp &matRcpp, int count, std::size_t phaseOne) {
    
    if (keepRes) {
        if (IsComb) {
            if (IsMultiset)
                MultisetComboResult(n, r, v, myReps, freqs, nRows, count, z, matRcpp, myFun);
            else
                ComboGenRes(n, r, v, IsRep, nRows, count, z, matRcpp, myFun);
        } else {
            if (IsMultiset)
                MultisetPermRes(n, r, v, nRows, count, z, matRcpp, myFun);
            else
                PermuteGenRes(n, r, v, IsRep, nRows, z, count, matRcpp, myFun);
        }
    } else {
        if (IsComb) {
            if (IsMultiset)
                MultisetCombination(n, r, v, myReps, freqs, count, nRows, z, matRcpp);
            else
                ComboGeneral(n, r, v, IsRep, count, nRows, z, matRcpp);
        } else {
            if (IsMultiset)
                MultisetPermutation(n, r, v, nRows, z, count, matRcpp);
            else if (generalReturn)
                PermuteGeneral(n, r, v, IsRep, nRows, z, count, matRcpp);
            else
                PermuteSerialDriver(n, r, v, IsRep, nRows, phaseOne, z, matRcpp);
        }
    }
}

// This is called when we can't easily produce a monotonic sequence overall,
// and we must generate and test every possible combination/permutation
template <typename typeRcpp, typename typeVector>
typeRcpp SpecCaseRet(int n, int m, std::vector<typeVector> v, bool IsRep, int nRows, 
                     bool keepRes, std::vector<int> z, double lower, std::string mainFun, 
                     bool IsMult, double computedRows, std::vector<std::string> compFunVec,
                     std::vector<typeVector> targetVals, bool IsComb, std::vector<int> myReps,
                     std::vector<int> freqs, bool bLower, double userRows) {
    
    if (!bLower) {
        if (computedRows > std::numeric_limits<int>::max())
            Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
        
        nRows = static_cast<int>(computedRows);
    }
    
    std::vector<typeVector> rowVec(m);
    std::vector<int> indexMatch;
    indexMatch.reserve(nRows);
    
    Rcpp::XPtr<funcPtr<typeVector>> xpFun = putFunPtrInXPtr<typeVector>(mainFun);
    funcPtr<typeVector> myFun = *xpFun;
    typeRcpp matRes = Rcpp::no_init_matrix(nRows, m + 1);
    
    // We pass keepRes = true (second true) as we need the results to determine which
    // results are within the constraint. The actual value of keepRes is utilized
    // below for the return matrix. The variable permNonTrivial, has no affect when
    // keepRes = true, so we pass it arbitrarily as true.
    GeneralReturn(n, m, v, IsRep, nRows, IsComb, myReps, 
                  freqs, z, true, IsMult, myFun, true, matRes, 0, 0);
    
    Rcpp::XPtr<compPtr<typeVector>> xpComp = putCompPtrInXPtr<typeVector>(compFunVec[0]);
    compPtr<typeVector> myComp = *xpComp;
    
    Rcpp::XPtr<compPtr<typeVector>> xpComp2 = xpComp;
    compPtr<typeVector> myComp2;
    
    if (compFunVec.size() == 1) {
        for (int i = 0; i < nRows; ++i) {
            const typeVector testVal = matRes(i, m);
            
            if (myComp(testVal, targetVals))
                indexMatch.push_back(i);
        }
    } else {
        xpComp2 = putCompPtrInXPtr<typeVector>(compFunVec[1]);
        myComp2 = *xpComp2;
        std::vector<typeVector> targetVals2 = targetVals;
        targetVals2.erase(targetVals2.begin());
        
        for (int i = 0; i < nRows; ++i) {
            const typeVector testVal = matRes(i, m);
            
            if (myComp(testVal, targetVals) || myComp2(testVal, targetVals2))
                indexMatch.push_back(i);
        }
    }
    
    const int numCols = (keepRes) ? (m + 1) : m;
    const int numMatches = indexMatch.size();
    
    if (bLower)
        nRows = numMatches;
    else
        nRows  = (numMatches > userRows) ? userRows : numMatches;
    
    typeRcpp returnMatrix = Rcpp::no_init_matrix(nRows, numCols);
    const int lastCol = keepRes ? (m + 1) : m;
    
    for (int i = 0; i < nRows; ++i)
        for (int j = 0; j < lastCol; ++j)
            returnMatrix(i, j) = matRes(indexMatch[i], j);
    
    return returnMatrix;
}

template <typename typeVector>
inline void BruteNextElem(int &ind, int lowBnd, typeVector targetMin, typeVector targetMax,
                          typeVector partial, std::size_t uR, const std::vector<typeVector> &v,
                          partialPtr<typeVector> partialFun) {
    
    typeVector dist = targetMax - partialFun(partial, v[ind], uR);
    int origInd = ind;
    
    while (ind > lowBnd && dist < 0) {
        --ind;
        dist = targetMax - partialFun(partial, v[ind], uR);
    }
    
    if ((targetMin - partialFun(partial, v[ind], uR)) > 0 && origInd != ind) {++ind;}
}

template <typename typeVector>
typeVector PartialReduce(std::size_t uR, typeVector partial, 
                         typeVector w, std::string myFun) {
    
    if (myFun == "prod") {
        partial /= w;
    } else if (myFun == "sum") {
        partial -= w;
    } else if (myFun == "mean") {
        partial = (partial * static_cast<double>(uR) - w) / static_cast<double>(uR - 1);
    }
    
    return partial;
}

template <typename typeVector>
int GetLowerBound(int n, int r, std::vector<typeVector> &v, bool isRep, bool isMult,
                  std::vector<int> &z, std::vector<int> &freqs, std::vector<typeVector> targetVals,
                  const std::vector<int> &Reps, funcPtr<typeVector> constraintFun,
                  partialPtr<typeVector> partialFun, std::string myFun) {
    
    const int lastElem = n - 1;
    const int lastCol = r - 1;
    const std::size_t uR = r;
    
    std::vector<typeVector> vPass(r);
    const typeVector targetMin = *std::min_element(targetVals.cbegin(), targetVals.cend());
    const typeVector targetMax = *std::max_element(targetVals.cbegin(), targetVals.cend());
    
    if (isRep) {
        std::fill(vPass.begin(), vPass.end(), v.back());
    } else if (isMult) {
        const int lenMinusR = freqs.size() - r;
        
        for (int i = freqs.size() - 1, j = 0; i >= lenMinusR; --i, ++j)
            vPass[j] = v[freqs[i]];
    } else {
        vPass.assign(v.crbegin(), v.crbegin() + r);
    }
    
    typeVector partial = constraintFun(vPass, uR - 1);
    const typeVector testMax = partialFun(partial, vPass.back(), uR);
    if (testMax < targetMin)  {return Partitions::noSoln;}

    if (isRep) {
        std::fill(vPass.begin(), vPass.end(), v[0]);
    } else if (isMult) {
        for (int i = 0; i < r; ++i)
            vPass[i] = v[freqs[i]];
    } else {
        vPass.assign(v.cbegin(), v.cbegin() + r);
    }

    const typeVector testMin = constraintFun(vPass, uR);
    if (testMin > targetMax)  {return Partitions::noSoln;}

    int zExpCurrPos = (isMult) ? freqs.size() - r : 0;
    int currPos = (isMult) ? freqs[zExpCurrPos] : ((isRep) ? lastElem : (n - r));
    int ind = currPos;

    int lowBnd = 0;
    std::vector<int> repsCounter;
    
    if (isMult)
        repsCounter.assign(Reps.cbegin(), Reps.cend());
    
    for (int i = 0; i < lastCol; ++i) {
        BruteNextElem(ind, lowBnd, targetMin, targetMax, partial, uR, v, partialFun);
        z[i] = ind;
        partial = partialFun(partial, v[ind], uR);
        
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
        ind = currPos;
        partial = PartialReduce(r, partial, v[currPos], myFun);
    }

    BruteNextElem(ind, lowBnd, targetMin, targetMax, partial, uR, v, partialFun);
    z[lastCol] = ind;
    return Partitions::solnExists;
}

// This function applys a constraint function to a vector v with respect
// to a constraint value "target". The main idea is that combinations are added
// successively, until a particular combination exceeds the given constraint
// value for a given constraint function. After this point, we can safely skip
// several combinations knowing that they will exceed the given constraint value.

template <typename typeRcpp, typename typeVector>
typeRcpp CombinatoricsConstraints(int n, int r, std::vector<typeVector> &v, bool isRep, 
                                  std::string myFun, std::vector<std::string> comparison,
                                  std::vector<typeVector> targetVals, double numRows, bool isComb,
                                  bool xtraCol, std::vector<int> &Reps, bool isMult, bool bUserRows) {
    
    // myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max";
    // The comparison vector is a comparison operator: 
    //                             "<", "<=", ">", ">=", "==", ">,<", ">=,<", ">,<=", ">=,<=";
    
    typeVector testVal;
    std::size_t count = 0;
    const bool partitionEsque = (comparison[0] == "==" && n > 1 && r > 1 && myFun != "max" && myFun != "min");
    const std::size_t maxRows = std::min(static_cast<double>(std::numeric_limits<int>::max()), numRows);
    
    const std::size_t uR = r;
    const std::size_t uR1 = r - 1;
    std::vector<typeVector> combinatoricsVec;
    std::vector<typeVector> resultsVec;
    
    if (bUserRows) {
        combinatoricsVec.reserve(uR * maxRows);
        resultsVec.reserve(maxRows);
    }
    
    Rcpp::XPtr<funcPtr<typeVector>> xpFun = putFunPtrInXPtr<typeVector>(myFun);
    funcPtr<typeVector> constraintFun = *xpFun;
    
    Rcpp::XPtr<partialPtr<typeVector>> xpPartial = putPartialPtrInXPtr<typeVector>(myFun);
    partialPtr<typeVector> partialFun = *xpPartial;
    
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
            
            auto itComp = std::find(compSpecial.cbegin(), compSpecial.cend(), comparison[nC]);
            
            if (itComp != compSpecial.end()) {
                int myIndex = std::distance(compSpecial.cbegin(), itComp);
                Rcpp::XPtr<compPtr<typeVector>> xpCompThree = putCompPtrInXPtr<typeVector>(compHelper[myIndex]);
                comparisonFunTwo = *xpCompThree;
            } else {
                comparisonFunTwo = *xpCompOne;
            }
        }
        
        std::vector<int> z(r);
        std::vector<typeVector> testVec(r);
        bool t_0 = true;
        bool t_1 = true;
        int myStart, maxZ = n - 1;
        const int r1 = r - 1;
        const int r2 = r - 2;
        
        if (r == 1) {
            int ind = 0;
            testVal = v[ind];
            t_0 = comparisonFunTwo(testVal, targetVals);
            
            while (t_0 && t_1) {
                if (comparisonFunOne(testVal, targetVals)) {
                    for (int k = 0; k < r; ++k)
                        combinatoricsVec.push_back(v[ind]);
                    
                    ++count;
                    
                    if (xtraCol)
                        resultsVec.push_back(testVal);
                    
                    t_1 =  (count < maxRows);
                }
                
                t_0 = ind != maxZ;
                
                if (t_0) {
                    ++ind;
                    testVal = v[ind];
                    t_0 = comparisonFunTwo(testVal, targetVals);
                }
            }
        } else {
        
            if (isMult) {
                int zExpSize = std::accumulate(Reps.cbegin(), Reps.cend(), 0);
                std::vector<int> zExpand, zIndex, zGroup(r), zPerm(r);
                
                for (int i = 0, k = 0; i < n; ++i) {
                    zIndex.push_back(k);
                    
                    for (int j = 0; j < Reps[i]; ++j, ++k)
                        zExpand.push_back(i);
                }
                
                if (partitionEsque) {
                    t_1 = GetLowerBound<typeVector>(n, r, v, isRep, isMult, z, zExpand, 
                                                    targetVals, Reps, constraintFun, partialFun, myFun);
                } else {
                    z.assign(zExpand.cbegin(), zExpand.cbegin() + r);
                }
                
                while (t_1) {
                    for (int i = 0; i < r; ++i)
                        testVec[i] = v[zExpand[zIndex[z[i]]]];
                    
                    const typeVector partialVal = constraintFun(testVec, uR1);
                    testVal = partialFun(partialVal, testVec.back(), uR);
                    t_0 = comparisonFunTwo(testVal, targetVals);
                    
                    while (t_0 && t_1) {
                        if (comparisonFunOne(testVal, targetVals)) {
                            myStart = count;
                            
                            if (isComb) {
                                for (int k = 0; k < r; ++k)
                                    combinatoricsVec.push_back(v[zExpand[zIndex[z[k]]]]);
                                
                                ++count;
                            } else {
                                for (int k = 0; k < r; ++k)
                                    zPerm[k] = zExpand[zIndex[z[k]]];
                                
                                do {
                                    for (int k = 0; k < r; ++k)
                                        combinatoricsVec.push_back(v[zPerm[k]]);
                                    
                                    ++count;
                                } while (std::next_permutation(zPerm.begin(), zPerm.end()) && count < maxRows);
                            }
                            
                            if (xtraCol)
                                for (std::size_t i = myStart; i < count; ++i)
                                    resultsVec.push_back(testVal);
                            
                            t_1 = count < maxRows;
                        }
                        
                        t_0 = z[r1] != maxZ;
                        
                        if (t_0) {
                            ++z[r1];
                            testVec[r1] = v[zExpand[zIndex[z[r1]]]];
                            testVal = partialFun(partialVal, testVec.back(), uR);
                            t_0 = comparisonFunTwo(testVal, targetVals);
                        }
                    }
                    
                    if (t_1) {
                        bool noChange = true;
                        
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
                                t_0 = comparisonFunTwo(testVal, targetVals);
                                noChange = false;
                                if (t_0) {break;}
                            }
                        }
                        
                        t_1 = (!noChange && t_0);
                    }
                }
                
            } else if (isRep) {
                
                v.erase(std::unique(v.begin(), v.end()), v.end());
                maxZ = static_cast<int>(v.size()) - 1;
                
                if (partitionEsque) {
                    std::vector<int> trivialFreqs;
                    t_1 = GetLowerBound<typeVector>(n, r, v, isRep, isMult, z, trivialFreqs, 
                                                    targetVals, Reps, constraintFun, partialFun, myFun);
                } else {
                    z.assign(r, 0);
                }
                
                while (t_1) {
                    for (int i = 0; i < r; ++i)
                        testVec[i] = v[z[i]];
                    
                    const typeVector partialVal = constraintFun(testVec, uR1);
                    testVal = partialFun(partialVal, testVec.back(), uR);
                    t_0 =comparisonFunTwo(testVal, targetVals);
                    
                    while (t_0 && t_1) {
                        if (comparisonFunOne(testVal, targetVals)) {
                            myStart = count;
                            
                            if (isComb) {
                                for (int k = 0; k < r; ++k)
                                    combinatoricsVec.push_back(v[z[k]]);
                                
                                ++count;
                            } else {
                                do {
                                    for (int k = 0; k < r; ++k)
                                        combinatoricsVec.push_back(v[z[k]]);
                                    
                                    ++count;
                                } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
                            }
                            
                            if (xtraCol)
                                for (std::size_t i = myStart; i < count; ++i)
                                    resultsVec.push_back(testVal);
                            
                            t_1 = count < maxRows;
                        }
                        
                        t_0 = z[r1] != maxZ;
                        
                        if (t_0) {
                            ++z[r1];
                            testVec[r1] = v[z[r1]];
                            testVal = partialFun(partialVal, testVec.back(), uR);
                            t_0 = comparisonFunTwo(testVal, targetVals);
                        }
                    }
                    
                    if (t_1) {
                        bool noChange = true;
                        
                        for (int i = r2; i >= 0; --i) {
                            if (z[i] != maxZ) {
                                ++z[i];
                                testVec[i] = v[z[i]];
                                
                                for (int k = (i+1); k < r; ++k) {
                                    z[k] = z[k-1];
                                    testVec[k] = v[z[k]];
                                }
                                
                                testVal = constraintFun(testVec, uR);
                                t_0 = comparisonFunTwo(testVal, targetVals);
                                noChange = false;
                                if (t_0) {break;}
                            }
                        }
                        
                        t_1 = (!noChange && t_0);
                    }
                }
                
            } else {
                
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
                
                if (partitionEsque) {
                    std::vector<int> trivialFreqs;
                    t_1 = GetLowerBound<typeVector>(n, r, v, isRep, isMult, z, trivialFreqs, 
                                                    targetVals, Reps, constraintFun, partialFun, myFun);
                } else {
                    std::iota(z.begin(), z.end(), 0);
                }
                
                while (t_1) {
                    for (int i = 0; i < r; ++i)
                        testVec[i] = v[z[i]];
                    
                    const typeVector partialVal = constraintFun(testVec, uR1);
                    testVal = partialFun(partialVal, testVec.back(), uR);
                    t_0 = comparisonFunTwo(testVal, targetVals);
                    
                    while (t_0 && t_1) {
                        if (comparisonFunOne(testVal, targetVals)) {
                            myStart = count;
                            
                            if (isComb) {
                                for (int k = 0; k < r; ++k)
                                    combinatoricsVec.push_back(v[z[k]]);
                                
                                ++count;
                            } else {
                                if (indexRows + count > maxRows)
                                    indexRows = maxRows - count;
                                    
                                for (int j = 0, myRow = 0; j < indexRows; ++j, myRow += r)
                                    for (int k = 0; k < r; ++k)
                                        combinatoricsVec.push_back(v[z[indexMatrix[myRow + k]]]);
                                
                                count += indexRows;
                            }
                            
                            if (xtraCol)
                                for (std::size_t i = myStart; i < count; ++i)
                                    resultsVec.push_back(testVal);
                            
                            t_1 = count < maxRows;
                        }
                        
                        t_0 = z[r1] != maxZ;
                        
                        if (t_0) {
                            ++z[r1];
                            testVec[r1] = v[z[r1]];
                            testVal = partialFun(partialVal, testVec.back(), uR);
                            t_0 = comparisonFunTwo(testVal, targetVals);
                        }
                    }
                    
                    if (t_1) {
                        bool noChange = true;
                        
                        for (int i = r2; i >= 0; --i) {
                            if (z[i] != (nMinusR + i)) {
                                ++z[i];
                                testVec[i] = v[z[i]];
                                
                                for (int k = (i + 1); k < r; ++k) {
                                    z[k] = z[k - 1] + 1;
                                    testVec[k] = v[z[k]];
                                }
                                
                                testVal = constraintFun(testVec, uR);
                                t_0 = comparisonFunTwo(testVal, targetVals);
                                noChange = false;
                                if (t_0) {break;}
                            }
                        }
                        
                        t_1 = (!noChange && t_0);
                    }
                }
            }
            
            targetVals.erase(targetVals.begin());
        }
    }
    const int numCols = xtraCol ? (r + 1) : r;
    typeRcpp combinatoricsMatrix = Rcpp::no_init_matrix(count, numCols);
    
    for (std::size_t i = 0, k = 0; i < count; ++i)
        for (int j = 0; j < r; ++j, ++k)
            combinatoricsMatrix(i, j) = combinatoricsVec[k];
    
    if (xtraCol)
        for (std::size_t i = 0; i < count; ++i)
            combinatoricsMatrix(i, r) = resultsVec[i];
    
    if (count > std::numeric_limits<int>::max()) {
        Rcpp::warning("The algorithm terminated early as the number of "
                          "results meeting the criteria exceeds 2^31 - 1.");
    }
    
    return combinatoricsMatrix;
}

#endif
