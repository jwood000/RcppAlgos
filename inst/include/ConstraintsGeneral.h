#ifndef CONSTRAINTS_GENERAL_H
#define CONSTRAINTS_GENERAL_H

#include "ConstraintsUtils.h"

// This function applys a constraint function to a vector v with respect
// to a constraint value "target". The main idea is that combinations are added
// successively, until a particular combination exceeds the given constraint
// value for a given constraint function. After this point, we can safely skip
// several combinations knowing that they will exceed the given constraint value.
template <typename typeRcpp, typename typeVector>
typeRcpp CombinatoricsConstraints(int n, int m, std::vector<typeVector> &v, bool isRep, 
                                  std::string myFun, std::vector<std::string> comparison,
                                  std::vector<typeVector> targetVals, double numRows, 
                                  bool isComb, bool xtraCol, std::vector<int> &Reps, 
                                  bool isMult, bool bUserRows, bool between) {
    
    // myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max";
    // The comparison vector is a comparison operator: 
    //                             "<", "<=", ">", ">=", "==", ">,<", ">=,<", ">,<=", ">=,<=";
    
    typeVector testVal;
    std::size_t count = 0;
    const bool partitionEsque = ((comparison[0] == "==" || between) && n > 1 && myFun != "max" && myFun != "min");
    const std::size_t maxRows = std::min(static_cast<double>(std::numeric_limits<int>::max()), numRows);
    
    std::vector<typeVector> combinatoricsVec;
    std::vector<typeVector> resultsVec;
    
    if (bUserRows) {
        combinatoricsVec.reserve(m * maxRows);
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
        
        std::vector<int> z(m);
        std::vector<typeVector> testVec(m);
        
        bool t_0 = true;
        bool t_1 = true;
        int myStart, maxZ = n - 1;
        
        const int m1 = m - 1;
        const int m2 = m - 2;
        
        if (m == 1) {
            int ind = 0;
            testVal = v[ind];
            t_0 = comparisonFunTwo(testVal, targetVals);
            
            while (t_0 && t_1) {
                if (comparisonFunOne(testVal, targetVals)) {
                    for (int k = 0; k < m; ++k)
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
                std::vector<int> zExpand, zIndex, zGroup(m), zPerm(m);
                
                for (int i = 0, k = 0; i < n; ++i) {
                    zIndex.push_back(k);
                    
                    for (int j = 0; j < Reps[i]; ++j, ++k)
                        zExpand.push_back(i);
                }
                
                if (partitionEsque) {
                    t_1 = GetLowerBound<typeVector>(n, m, v, isRep, isMult, z, zExpand, 
                                                    targetVals, Reps, constraintFun, partialFun, myFun);
                } else {
                    z.assign(zExpand.cbegin(), zExpand.cbegin() + m);
                }
                
                while (t_1) {
                    for (int i = 0; i < m; ++i)
                        testVec[i] = v[zExpand[zIndex[z[i]]]];
                    
                    const typeVector partialVal = constraintFun(testVec, m1);
                    testVal = partialFun(partialVal, testVec.back(), m);
                    t_0 = comparisonFunTwo(testVal, targetVals);
                    
                    while (t_0 && t_1) {
                        if (comparisonFunOne(testVal, targetVals)) {
                            myStart = count;
                            
                            if (isComb) {
                                for (int k = 0; k < m; ++k)
                                    combinatoricsVec.push_back(v[zExpand[zIndex[z[k]]]]);
                                
                                ++count;
                            } else {
                                for (int k = 0; k < m; ++k)
                                    zPerm[k] = zExpand[zIndex[z[k]]];
                                
                                do {
                                    for (int k = 0; k < m; ++k)
                                        combinatoricsVec.push_back(v[zPerm[k]]);
                                    
                                    ++count;
                                } while (std::next_permutation(zPerm.begin(), zPerm.end()) && count < maxRows);
                            }
                            
                            if (xtraCol)
                                for (std::size_t i = myStart; i < count; ++i)
                                    resultsVec.push_back(testVal);
                            
                            t_1 = count < maxRows;
                        }
                        
                        t_0 = z[m1] != maxZ;
                        
                        if (t_0) {
                            ++z[m1];
                            testVec[m1] = v[zExpand[zIndex[z[m1]]]];
                            testVal = partialFun(partialVal, testVec.back(), m);
                            t_0 = comparisonFunTwo(testVal, targetVals);
                        }
                    }
                    
                    if (t_1) {
                        bool noChange = true;
                        
                        for (int i = m2; i >= 0; --i) {
                            if (zExpand[zIndex[z[i]]] != zExpand[zExpSize - m + i]) {
                                ++z[i];
                                testVec[i] = v[zExpand[zIndex[z[i]]]];
                                zGroup[i] = zIndex[z[i]];
                                
                                for (int k = (i+1); k < m; ++k) {
                                    zGroup[k] = zGroup[k-1] + 1;
                                    z[k] = zExpand[zGroup[k]];
                                    testVec[k] = v[zExpand[zIndex[z[k]]]];
                                }
                                
                                testVal = constraintFun(testVec, m);
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
                    t_1 = GetLowerBound<typeVector>(n, m, v, isRep, isMult, z, trivialFreqs, 
                                                    targetVals, Reps, constraintFun, partialFun, myFun);
                } else {
                    z.assign(m, 0);
                }
                
                while (t_1) {
                    for (int i = 0; i < m; ++i)
                        testVec[i] = v[z[i]];
                    
                    const typeVector partialVal = constraintFun(testVec, m1);
                    testVal = partialFun(partialVal, testVec.back(), m);
                    t_0 =comparisonFunTwo(testVal, targetVals);
                    
                    while (t_0 && t_1) {
                        if (comparisonFunOne(testVal, targetVals)) {
                            myStart = count;
                            
                            if (isComb) {
                                for (int k = 0; k < m; ++k)
                                    combinatoricsVec.push_back(v[z[k]]);
                                
                                ++count;
                            } else {
                                do {
                                    for (int k = 0; k < m; ++k)
                                        combinatoricsVec.push_back(v[z[k]]);
                                    
                                    ++count;
                                } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
                            }
                            
                            if (xtraCol)
                                for (std::size_t i = myStart; i < count; ++i)
                                    resultsVec.push_back(testVal);
                            
                            t_1 = count < maxRows;
                        }
                        
                        t_0 = z[m1] != maxZ;
                        
                        if (t_0) {
                            ++z[m1];
                            testVec[m1] = v[z[m1]];
                            testVal = partialFun(partialVal, testVec.back(), m);
                            t_0 = comparisonFunTwo(testVal, targetVals);
                        }
                    }
                    
                    if (t_1) {
                        bool noChange = true;
                        
                        for (int i = m2; i >= 0; --i) {
                            if (z[i] != maxZ) {
                                ++z[i];
                                testVec[i] = v[z[i]];
                                
                                for (int k = (i+1); k < m; ++k) {
                                    z[k] = z[k-1];
                                    testVec[k] = v[z[k]];
                                }
                                
                                testVal = constraintFun(testVec, m);
                                t_0 = comparisonFunTwo(testVal, targetVals);
                                noChange = false;
                                if (t_0) {break;}
                            }
                        }
                        
                        t_1 = (!noChange && t_0);
                    }
                }
                
            } else {
                
                const int nMinusM = (n - m);
                int indexRows = isComb ? 0 : static_cast<int>(NumPermsNoRep(m, m1));
                auto indexMatrix = FromCpp14::make_unique<int[]>(indexRows * m);
                
                if (!isComb) {
                    indexRows = static_cast<int>(NumPermsNoRep(m, m1));
                    std::vector<int> indexVec(m);
                    std::iota(indexVec.begin(), indexVec.end(), 0);
                    
                    for (int i = 0, myRow = 0; i < indexRows; ++i, myRow += m) {
                        for (int j = 0; j < m; ++j)
                            indexMatrix[myRow + j] = indexVec[j];
                        
                        std::next_permutation(indexVec.begin(), indexVec.end());
                    }
                }
                
                if (partitionEsque) {
                    std::vector<int> trivialFreqs;
                    t_1 = GetLowerBound<typeVector>(n, m, v, isRep, isMult, z, trivialFreqs, 
                                                    targetVals, Reps, constraintFun, partialFun, myFun);
                } else {
                    std::iota(z.begin(), z.end(), 0);
                }
                
                while (t_1) {
                    for (int i = 0; i < m; ++i)
                        testVec[i] = v[z[i]];
                    
                    const typeVector partialVal = constraintFun(testVec, m1);
                    testVal = partialFun(partialVal, testVec.back(), m);
                    t_0 = comparisonFunTwo(testVal, targetVals);
                    
                    while (t_0 && t_1) {
                        if (comparisonFunOne(testVal, targetVals)) {
                            myStart = count;
                            
                            if (isComb) {
                                for (int k = 0; k < m; ++k)
                                    combinatoricsVec.push_back(v[z[k]]);
                                
                                ++count;
                            } else {
                                if (indexRows + count > maxRows)
                                    indexRows = maxRows - count;
                                    
                                for (int j = 0, myRow = 0; j < indexRows; ++j, myRow += m)
                                    for (int k = 0; k < m; ++k)
                                        combinatoricsVec.push_back(v[z[indexMatrix[myRow + k]]]);
                                
                                count += indexRows;
                            }
                            
                            if (xtraCol)
                                for (std::size_t i = myStart; i < count; ++i)
                                    resultsVec.push_back(testVal);
                            
                            t_1 = count < maxRows;
                        }
                        
                        t_0 = z[m1] != maxZ;
                        
                        if (t_0) {
                            ++z[m1];
                            testVec[m1] = v[z[m1]];
                            testVal = partialFun(partialVal, testVec.back(), m);
                            t_0 = comparisonFunTwo(testVal, targetVals);
                        }
                    }
                    
                    if (t_1) {
                        bool noChange = true;
                        
                        for (int i = m2; i >= 0; --i) {
                            if (z[i] != (nMinusM + i)) {
                                ++z[i];
                                testVec[i] = v[z[i]];
                                
                                for (int k = (i + 1); k < m; ++k) {
                                    z[k] = z[k - 1] + 1;
                                    testVec[k] = v[z[k]];
                                }
                                
                                testVal = constraintFun(testVec, m);
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
    const int numCols = xtraCol ? (m + 1) : m;
    typeRcpp combinatoricsMatrix = Rcpp::no_init_matrix(count, numCols);
    
    for (std::size_t i = 0, k = 0; i < count; ++i)
        for (int j = 0; j < m; ++j, ++k)
            combinatoricsMatrix(i, j) = combinatoricsVec[k];
    
    if (xtraCol)
        for (std::size_t i = 0; i < count; ++i)
            combinatoricsMatrix(i, m) = resultsVec[i];
    
    if (count > std::numeric_limits<int>::max()) {
        Rcpp::warning("The algorithm terminated early as the number of "
                          "results meeting the criteria exceeds 2^31 - 1.");
    }
    
    return combinatoricsMatrix;
}

#endif
