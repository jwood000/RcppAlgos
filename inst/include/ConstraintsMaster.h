#ifndef CONSTRAINTS_MASTER_H
#define CONSTRAINTS_MASTER_H

#include "Combinations.h"
#include "Permutations.h"
#include "ComboResults.h"
#include "PermuteResults.h"

const std::vector<std::string> compForms = {"<", ">", "<=", ">=", "==", "=<", "=>"};
const std::vector<std::string> compSpecial = {"==", ">,<", ">=,<", ">,<=", ">=,<="};
const std::vector<std::string> compHelper = {"<=", "<", "<", "<=", "<="};

template <typename typeRcpp, typename typeVector>
void GeneralReturn(int n, int m, std::vector<typeVector> v, bool IsRep, int nRows, bool IsComb,
                   std::vector<int> myReps, std::vector<int> freqs, std::vector<int> z,
                   bool permNonTriv, bool IsMultiset, funcPtr<typeVector> myFun,
                   bool keepRes, typeRcpp &matRcpp, int count) {
    if (keepRes) {
        if (IsComb) {
            if (IsMultiset)
                MultisetComboResult(n, m, v, myReps, freqs, nRows, count, z, matRcpp, myFun);
            else
                ComboGenRes(n, m, v, IsRep, nRows, count, z, matRcpp, myFun);
        } else {
            if (IsMultiset)
                MultisetPermRes(n, m, v, nRows, count, z, matRcpp, myFun);
            else
                PermuteGenRes(n, m, v, IsRep, nRows, z, count, permNonTriv, matRcpp, myFun);
        }
    } else {
        if (IsComb) {
            if (IsMultiset)
                MultisetCombination(n, m, v, myReps, freqs, count, nRows, z, matRcpp);
            else
                ComboGeneral(n, m, v, IsRep, count, nRows, z, matRcpp);
        } else {
            if (IsMultiset)
                MultisetPermutation(n, m, v, nRows, z, count, matRcpp);
            else
                PermuteGeneral(n, m, v, IsRep, nRows, z, count, permNonTriv, matRcpp);
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
                     std::vector<int> freqs, bool bLower, bool permNT, double userRows, double tol) {
    
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
    typeVector testVal;
    
    GeneralReturn(n, m, v, IsRep, nRows, IsComb, myReps, 
                  freqs, z, permNT, IsMult, myFun, true, matRes, 0);
    
    Rcpp::XPtr<compPtr<typeVector>> xpComp = putCompPtrInXPtr<typeVector>(compFunVec[0]);
    compPtr<typeVector> myComp = *xpComp;
    
    Rcpp::XPtr<compPtr<typeVector>> xpComp2 = xpComp;
    compPtr<typeVector> myComp2;
    
    if (compFunVec.size() == 1) {
        for (int i = 0; i < nRows; ++i) {
            testVal = matRes(i, m);
            
            if (myComp(testVal, targetVals))
                indexMatch.push_back(i);
        }
    } else {
        xpComp2 = putCompPtrInXPtr<typeVector>(compFunVec[1]);
        myComp2 = *xpComp2;
        std::vector<typeVector> targetVals2 = targetVals;
        targetVals2.erase(targetVals2.begin());
        
        for (int i = 0; i < nRows; ++i) {
            testVal = matRes(i, m);
            bool Success = myComp(testVal, targetVals) || myComp2(testVal, targetVals2);
            
            if (Success)
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

// This function applys a constraint function to a vector v with respect
// to a constraint value "target". The main idea is that combinations are added
// successively, until a particular combination exceeds the given constraint
// value for a given constraint function. After this point, we can safely skip
// several combinations knowing that they will exceed the given constraint value.

template <typename typeRcpp, typename typeVector>
typeRcpp CombinatoricsConstraints(int n, int r, std::vector<typeVector> &v, bool isRep, std::string myFun,
                                  std::vector<std::string> comparison, std::vector<typeVector> targetVals, double numRows,
                                  bool isComb, bool xtraCol, std::vector<int> &Reps, bool isMult, double tol, bool bUserRows) {
    
    // myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max";
    // The comparison vector is a comparison operator: 
    //                             "<", "<=", ">", ">=", "==", ">,<", ">=,<", ">,<=", ">=,<=";
    
    typeVector testVal;
    std::size_t count = 0;
    const std::size_t maxRows = std::min(static_cast<double>(std::numeric_limits<int>::max()), numRows);
    
    const std::size_t uR = r;
    std::vector<typeVector> combinatoricsVec;
    std::vector<typeVector> resultsVec;
    
    if (bUserRows) {
        combinatoricsVec.reserve(uR * maxRows);
        resultsVec.reserve(maxRows);
    }
    
    Rcpp::XPtr<funcPtr<typeVector>> xpFun = putFunPtrInXPtr<typeVector>(myFun);
    funcPtr<typeVector> constraintFun = *xpFun;
    
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
            
            z.assign(zExpand.cbegin(), zExpand.cbegin() + r);
            
            while (keepGoing) {
                
                t_2 = true;
                for (int i = 0; i < r; ++i)
                    testVec[i] = v[zExpand[zIndex[z[i]]]];
                
                testVal = constraintFun(testVec, uR);
                t = comparisonFunTwo(testVal, targetVals);
                
                while (t && t_2 && keepGoing) {
                    
                    testVal = constraintFun(testVec, uR);
                    t_1 = comparisonFunOne(testVal, targetVals);
                    
                    if (t_1) {
                        myStart = count;
                        
                        if (isComb) {
                            for (int k = 0; k < r; ++k)
                                combinatoricsVec.push_back(v[zExpand[zIndex[z[k]]]]);
                            
                            ++count;
                        } else {
                            for (int k = 0; k < r; ++k)
                                zPerm[k] = zExpand[zIndex[z[k]]];
                            
                            numIter = static_cast<int>(NumPermsWithRep(zPerm));
                            
                            if ((numIter + count) > maxRows)
                                numIter = maxRows - count;
                            
                            for (int i = 0; i < numIter; ++i) {
                                for (int k = 0; k < r; ++k)
                                    combinatoricsVec.push_back(v[zPerm[k]]);
                                
                                std::next_permutation(zPerm.begin(), zPerm.end());
                            }
                            
                            count += numIter;
                        }
                        
                        if (xtraCol)
                            for (std::size_t i = myStart; i < count; ++i)
                                resultsVec.push_back(testVal);
                    }
                    
                    keepGoing = (count < maxRows);
                    t_2 = (z[r1] != maxZ);
                    
                    if (t_2) {
                        ++z[r1];
                        testVec[r1] = v[zExpand[zIndex[z[r1]]]];
                        testVal = constraintFun(testVec, uR);
                        t = comparisonFunTwo(testVal, targetVals);
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
                            t = comparisonFunTwo(testVal, targetVals);
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
                t = comparisonFunTwo(testVal, targetVals);
                
                while (t && t_2 && keepGoing) {
                    
                    testVal = constraintFun(testVec, uR);
                    t_1 = comparisonFunOne(testVal, targetVals);
                    
                    if (t_1) {
                        myStart = count;
                        
                        if (isComb) {
                            for (int k = 0; k < r; ++k)
                                combinatoricsVec.push_back(v[z[k]]);
                            
                            ++count;
                        } else {
                            numIter = static_cast<int>(NumPermsWithRep(z));
                            
                            if ((numIter + count) > maxRows)
                                numIter = maxRows - count;
                            
                            for (int i = 0; i < numIter; ++i) {
                                for (int k = 0; k < r; ++k)
                                    combinatoricsVec.push_back(v[z[k]]);
                                
                                std::next_permutation(z.begin(), z.end());
                            }
                            
                            count += numIter;
                        }
                        
                        if (xtraCol)
                            for (std::size_t i = myStart; i < count; ++i)
                                resultsVec.push_back(testVal);
                        
                        keepGoing = (count < maxRows);
                    }
                    
                    t_2 = (z[r1] != maxZ);
                    
                    if (t_2) {
                        ++z[r1];
                        testVec[r1] = v[z[r1]];
                        testVal = constraintFun(testVec, uR);
                        t = comparisonFunTwo(testVal, targetVals);
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
                            t = comparisonFunTwo(testVal, targetVals);
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
                t = comparisonFunTwo(testVal, targetVals);
                
                while (t && t_2 && keepGoing) {
                    
                    testVal = constraintFun(testVec, uR);
                    t_1 = comparisonFunOne(testVal, targetVals);
                    
                    if (t_1) {
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
                        
                        keepGoing = (count < maxRows);
                    }
                    
                    t_2 = (z[r1] != maxZ);
                    
                    if (t_2) {
                        ++z[r1];
                        testVec[r1] = v[z[r1]];
                        testVal = constraintFun(testVec, uR);
                        t = comparisonFunTwo(testVal, targetVals);
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
                            t = comparisonFunTwo(testVal, targetVals);
                            if (t) {break;}
                        }
                    }
                    
                    if (!t || zCheck == z) {keepGoing = false;}
                }
            }
        }
        
        targetVals.erase(targetVals.begin());
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
