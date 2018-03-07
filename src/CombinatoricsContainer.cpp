#include <RcppAlgos.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
using namespace Rcpp;

template <typename typeRcpp>
typeRcpp SubMat(typeRcpp m, int n) {
    int k = m.ncol();
    typeRcpp subMatrix(n,k);
    for (int i = 0; i < n; i++) {subMatrix(i,_) = m(i,_);}
    return subMatrix;
}

// Comparison functions
bool lessCpp(double& x, double& y) {return x < y;}
bool greaterCpp(double& x, double& y) {return x > y;}
bool lessEqualCpp(double& x, double& y) {return x <= y;}
bool greaterEqualCpp(double& x, double& y) {return x >= y;}
bool equalCpp(double& x, double& y) {return x == y;}

typedef double (*funcPtr)(std::vector<double>& x);
typedef bool (*compPtr)(double& x, double& y);

XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
    if (fstr == "prod")
        return(XPtr<funcPtr>(new funcPtr(&prodCpp)));
    else if (fstr == "sum")
        return(XPtr<funcPtr>(new funcPtr(&sumCpp)));
    else if (fstr == "mean")
        return(XPtr<funcPtr>(new funcPtr(&meanCpp)));
    else if (fstr == "max")
        return(XPtr<funcPtr>(new funcPtr(&maxCpp)));
    else if (fstr == "min")
        return(XPtr<funcPtr>(new funcPtr(&minCpp)));
    else
        return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}

XPtr<compPtr> putCompPtrInXPtr(std::string fstr) {
    if (fstr == "<")
        return(XPtr<compPtr>(new compPtr(&lessCpp)));
    else if (fstr == ">")
        return(XPtr<compPtr>(new compPtr(&greaterCpp)));
    else if (fstr == "<=")
        return(XPtr<compPtr>(new compPtr(&lessEqualCpp)));
    else if (fstr == ">=")
        return(XPtr<compPtr>(new compPtr(&greaterEqualCpp)));
    else if (fstr == "==")
        return(XPtr<compPtr>(new compPtr(&equalCpp)));
    else
        return XPtr<compPtr>(R_NilValue); // runtime error as NULL no XPtr
}

// This function applys a constraint function to a vector v with respect
// to a constraint value "lim". The main idea is that combinations are added
// successively, until a particular combination exceeds the given constraint
// value for a given constraint function. After this point, we can safely skip
// several combinations knowing that they will exceed the given constraint value.
template <typename typeRcpp>
typeRcpp CombinatoricsConstraints(int n, int r, std::vector<double> v,
                                       bool repetition, std::string myFun, 
                                       std::string myComparison, double lim,
                                       int numRows, bool isComb, bool xtraCol,
                                       std::vector<int> Reps, bool isMult) {
    
    // where myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max"
    // and myComparison is a comparison operator: "<", "<=", ">", or ">="
    
    double testVal;
    int count = 0, vSize = v.size(), numCols;
    int testRows, i, j, numRows2 = 0;
    
    // constraintFun is a pointer to one of the functions defined above
    XPtr<funcPtr> xpFun = putFunPtrInXPtr(myFun);
    funcPtr constraintFun = *xpFun;
    
    // comparisonFun is a pointer to one of the comparisons defined above
    XPtr<compPtr> xpCompOne = putCompPtrInXPtr(myComparison);
    compPtr comparisonFunOne = *xpCompOne;
    
    XPtr<compPtr> xpCompTwo = xpCompOne;
    compPtr comparisonFunTwo;

    if (myComparison == ">" || myComparison == ">=") {
        if (isMult) {
            for (i = 0; i < (n-1); i++) {
                for (j = (i+1); j < n; j++) {
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
            for (i = 0; i < (n-1); i++) {
                for (j = (i+1); j < n; j++) {
                    if (v[i] > v[j]) {
                        std::swap(v[i], v[j]);
                        std::swap(Reps[i], Reps[j]);
                    }
                }
            }
        } else {
            std::sort(v.begin(), v.end());
        }
        
        if (myComparison == "==") {
            XPtr<compPtr> xpCompThree = putCompPtrInXPtr("<=");
            comparisonFunTwo = *xpCompThree;
        } else {
            comparisonFunTwo = *xpCompOne;
        }
    }
    
    // We have to take special care for the case when we have
    // permutations of a multiset. In the parent function
    // CombinatoricsRcpp, when this case occurs, we set
    // the number of rows equal to the number of rows (numRows) we would
    // need as if we were finding all COMBINATIONS of a multiset.
    // This is because we initially need to find all combinations
    // that meet the constraints and subsequently find the
    // permutations using std::next_permutation. There is logic
    // in the parent function that accounts for the user passing
    // a rowCap that exceeds numRows set above, as there are
    // more possible permutations than combinations. The following
    // captures the possible rowCap and correctly sets numRows to
    // the number of combinations of a multiset. We then use
    // these two numbers to determine when we need to terminate.
    
    bool executeCap = true;
    testRows = numRows;
    
    if (!isComb && isMult) {
        numRows = MultisetCombRowNum(n, r, Reps);
        if (testRows == numRows) {executeCap = false;}
    }
    
    if (xtraCol) {numCols = r+1;} else {numCols = r;}
    typeRcpp combinatoricsMatrix(numRows, numCols);
    
    std::vector<int> z, zCheck, zPerm(r);
    std::vector<double> testVec(r);
    bool t_1, t_2, t = true, keepGoing = true;
    int r1 = r - 1, r2 = r - 2, k = 0; 
    int numIter, maxZ = n - 1;
    
    if (isMult) {
        int zExpSize = 0;
        std::vector<int> eachRowCount;
        if (!isComb) {eachRowCount.reserve(numRows);}
        std::vector<int> zExpand, zIndex, zGroup(r);
        for (i = 0; i < n; i++) {zExpSize += Reps[i];}
        
        zIndex.reserve(n);
        zExpand.reserve(zExpSize);
        for (i = 0; i < n; i++) {
            zIndex.push_back(k);
            for (j = 0; j < Reps[i]; j++, k++) {
                zExpand.push_back(i);
            }
        }
        
        z.reserve(r);
        for (i = 0; i < r; i++) {
            z.push_back(zExpand[i]);
        }
        
        while (keepGoing) {
            
            t_2 = true;
            for (i = 0; i < r; i++) {testVec[i] = v[zExpand[zIndex[z[i]]]];}
            testVal = constraintFun(testVec);
            t = comparisonFunTwo(testVal, lim);
            
            while (t && t_2 && keepGoing) {
                
                testVal = constraintFun(testVec);
                t_1 = comparisonFunOne(testVal, lim);
                
                if (t_1) {
                    for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = v[zExpand[zIndex[z[k]]]];}
                    count++;
                    
                    if (!isComb) {
                        for (k=0; k < r; k++) {zPerm[k] = zExpand[zIndex[z[k]]];}
                        numIter = (int)NumPermsWithRep(zPerm);
                        numRows2 += numIter;
                        eachRowCount.push_back(numIter);
                    }
                    
                    if (executeCap) {
                        keepGoing = (count < testRows);    
                    }
                }
                
                t_2 = (z[r1] != maxZ);
                
                if (t_2) {
                    z[r1]++;
                    testVec[r1] = v[zExpand[zIndex[z[r1]]]];
                    testVal = constraintFun(testVec);
                    t = comparisonFunTwo(testVal, lim);
                }
            }
            
            if (keepGoing) {
                zCheck = z;
                for (i = r2; i >= 0; i--) {
                    if (zExpand[zIndex[z[i]]] != zExpand[zExpSize - r + i]) {
                        z[i]++;
                        testVec[i] = v[zExpand[zIndex[z[i]]]];
                        zGroup[i] = zIndex[z[i]];
                        for (k = (i+1); k < r; k++) {
                            zGroup[k] = zGroup[k-1] + 1;
                            z[k] = zExpand[zGroup[k]];
                            testVec[k] = v[zExpand[zIndex[z[k]]]];
                        }
                        testVal = constraintFun(testVec);
                        t = comparisonFunTwo(testVal, lim);
                        if (t) {break;}
                    }
                }
                if (!t) {
                    keepGoing = false;
                } else if (zCheck == z) {
                    keepGoing = false;
                }
            }
        }
  
        if (!isComb) {
            typeRcpp permuteMatrix(numRows2, numCols);
            std::vector<double>::iterator myIter, permEnd;
            int count2 = 0, segment = 0;
            bool bTooMany = false;
            
            for (i = 0; i < count; i++) {
                // populate permVec with next combination
                for (j = 0; j < r; j++) {testVec[j] = combinatoricsMatrix(i, j);}
                
                segment = eachRowCount[i];
                if (executeCap) {
                    if (eachRowCount[i] + count2 > testRows) {
                        segment = testRows - count2;
                        bTooMany = true;
                    }
                }
                
                // populate permuteMatrix with all permutations
                // of a particular combination
                for (j = 0; j < segment; j++, count2++) {
                    permEnd = testVec.end();
                    for (myIter = testVec.begin(), k = 0; myIter < permEnd; myIter++, k++) {
                        permuteMatrix(count2, k) = *myIter;
                    }
                    std::next_permutation(testVec.begin(), permEnd);
                }
                if (bTooMany) {break;}
            }
            return SubMat(permuteMatrix, count2);
        }
        
    } else if (repetition) {
        v.erase(std::unique(v.begin(), v.end()), v.end());
        vSize = v.size();
        z.assign(r, 0);
        maxZ = vSize - 1;
        
        while (keepGoing) {
            
            t_2 = true;
            for (i = 0; i < r; i++) {testVec[i] = v[z[i]];}
            testVal = constraintFun(testVec);
            t = comparisonFunTwo(testVal, lim);
            
            while (t && t_2 && keepGoing) {
                
                testVal = constraintFun(testVec);
                t_1 = comparisonFunOne(testVal, lim);
                
                if (t_1) {
                    if (isComb) {
                        for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = v[z[k]];}
                        count++;
                    } else {
                        zPerm = z;
                        numIter = (int)NumPermsWithRep(zPerm);
                        if (numIter + count > numRows) {numIter = numRows - count;}
                        for (i=0; i < numIter; i++) {
                            for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = v[zPerm[k]];}
                            std::next_permutation(zPerm.begin(), zPerm.end());
                            count++;
                        }
                    }
                    keepGoing = (count < numRows);
                }
                
                t_2 = (z[r1] != maxZ);
                
                if (t_2) {
                    z[r1]++;
                    testVec[r1] = v[z[r1]];
                    testVal = constraintFun(testVec);
                    t = comparisonFunTwo(testVal, lim);
                }
            }
            
            if (keepGoing) {
                zCheck = z;
                for (i = r2; i >= 0; i--) {
                    if (z[i] != maxZ) {
                        z[i]++;
                        testVec[i] = v[z[i]];
                        for (k = (i+1); k < r; k++) {
                            z[k] = z[k-1];
                            testVec[k] = v[z[k]];
                        }
                        testVal = constraintFun(testVec);
                        t = comparisonFunTwo(testVal, lim);
                        if (t) {break;}
                    }
                }
                if (!t) {
                    keepGoing = false;
                } else if (zCheck == z) {
                    keepGoing = false;
                }
            }
        }
    } else {
        
        for (i = 0; i < r; i++) {z.push_back(i);}
        int indexRows = 0, nMinusR = (n - r);
        IntegerMatrix indexMatrix;
        if (!isComb) {
            indexRows = (int)NumPermsNoRep(r, r1);
            indexMatrix = MakeIndexHeaps(indexRows, r);
        }

        while (keepGoing) {
            
            t_2 = true;
            for (i = 0; i < r; i++) {testVec[i] = v[z[i]];}
            testVal = constraintFun(testVec);
            t = comparisonFunTwo(testVal, lim);
            
            while (t && t_2 && keepGoing) {
                
                testVal = constraintFun(testVec);
                t_1 = comparisonFunOne(testVal, lim);
                
                if (t_1) {
                    if (isComb) {
                        for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = v[z[k]];}
                        count++;
                    } else {
                        if (indexRows + count > numRows) {indexRows = numRows - count;}
                        for (j = 0; j < indexRows; j++) {
                            for (k = 0; k < r; k++) {
                                combinatoricsMatrix(count, k) = v[z[indexMatrix(j, k)]];
                            }
                            count++;
                        }
                    }
                    keepGoing = (count < numRows);
                }
                
                t_2 = (z[r1] != maxZ);
                
                if (t_2) {
                    z[r1]++;
                    testVec[r1] = v[z[r1]];
                    testVal = constraintFun(testVec);
                    t = comparisonFunTwo(testVal, lim);
                }
            }
            
            if (keepGoing) {
                zCheck = z;
                for (i = r2; i >= 0; i--) {
                    if (z[i] != (nMinusR + i)) {
                        z[i]++;
                        testVec[i] = v[z[i]];
                        for (k = (i+1); k < r; k++) {
                            z[k] = z[k-1]+1;
                            testVec[k] = v[z[k]];
                        }
                        testVal = constraintFun(testVec);
                        t = comparisonFunTwo(testVal, lim);
                        if (t) {break;}
                    }
                }
                if (!t) {
                    keepGoing = false;
                } else if (zCheck == z) {
                    keepGoing = false;
                }
            }
        }
    }
        
    return(SubMat(combinatoricsMatrix, count));
}

template <typename typeRcpp, typename typeStd>
typeRcpp ComboGeneral(int n, int r, std::vector<typeStd> v,
                      bool repetition, int numRows, bool xtraCol) {
    std::sort(v.begin(), v.end());
    std::vector<int> z;
    int r1 = r - 1, r2 = r - 2;
    int k, i, numIter, vSize;
    int numCols, maxZ, count = 0;
    bool needsSubsetting = false;
    if (xtraCol) {numCols  = r + 1;} else {numCols = r;}
    typeRcpp combinationMatrix(numRows, numCols);
    
    if (repetition) {
        v.erase(std::unique(v.begin(), v.end()), v.end());
        vSize = v.size();
        int testRows = GetRowNum(vSize, r);
        if (testRows < numRows) {
            needsSubsetting = true;
            numRows = testRows;
        }
        z.assign(r, 0);
        maxZ = vSize - 1;

        while (count < numRows) {
            numIter = vSize - z[r1];
            if (numIter + count > numRows) {numIter = numRows - count;}
            for (i = 0; i < numIter; i++) {
                for (k=0; k < r; k++) {combinationMatrix(count, k) = v[z[k]];}
                count++;
                z[r1]++;
            }

            for (i = r2; i >= 0; i--) {
                if (z[i] != maxZ) {
                    z[i]++;
                    for (k = (i+1); k < r; k++) {z[k] = z[k-1];}
                    break;
                }
            }
        }
    } else {
        for (i = 0; i < r; i++) {z.push_back(i);}
        while (count < numRows) {
            numIter = n - z[r1];
            if (numIter + count > numRows) {numIter = numRows - count;}
            for (i = 0; i < numIter; i++) {
                for (k = 0; k < r; k++) {combinationMatrix(count, k) = v[z[k]];}
                count++;
                z[r1]++;
            }

            for (i = r2; i >= 0; i--) {
                if (z[i] != (n - r + i)) {
                    z[i]++;
                    for (k = (i+1); k < r; k++) {z[k] = z[k-1]+1;}
                    break;
                }
            }
        }
    }
    
    if (needsSubsetting) {
        return SubMat(combinationMatrix, count);
    } else {
        return combinationMatrix;
    }
}

template <typename typeRcpp, typename typeStd>
typeRcpp PermuteGeneral(int n, int r, std::vector<typeStd> v,
                        bool repetition, int numRows, bool xtraCol) {
    unsigned long int uN = n, uR = r, uRowN = numRows;
    unsigned long int i = 0, j, k, chunk, numCols, segment;
    if (xtraCol) {numCols  = uR + 1;} else {numCols = uR;}
    typeRcpp permuteMatrix(uRowN, numCols);
    bool bTooMany = false;
    
    if (repetition) {
        unsigned long int groupLen = 1, repLen = 1, myCol;
        typename std::vector<typeStd>::iterator m, vBeg, vEnd;
        vBeg = v.begin(); vEnd = v.end();
        for (i = uR; i > 0; i--) {
            myCol = i-1;
            groupLen *= uN;
            chunk = 0;
            bTooMany = false;
            segment = repLen;
            for (k = 0; k < uRowN; k += groupLen) {
                for (m = vBeg; m < vEnd; m++) {
                    if (chunk + repLen > uRowN) {
                        segment = numRows - chunk;
                        bTooMany = true;
                    }
                    for (j = 0; j < segment; j++) {permuteMatrix(chunk + j, myCol) = *m;}
                    if (bTooMany) {break;}
                    chunk += repLen;
                }
                if (bTooMany) {break;}
            }
            repLen *= uN;
        }
    } else {
        unsigned long int combRows = (unsigned long int)nChooseK(n, r);
        typeRcpp myCombs = ComboGeneral<typeRcpp>(n,r,v,false,combRows,false);
        unsigned long int indexRows = (unsigned long int)NumPermsNoRep(r, r-1);
        IntegerMatrix indexMatrix = MakeIndexHeaps(indexRows, uR);

        chunk = 0;
        for (i = 0; i < combRows; i++) {
            if (chunk + indexRows > uRowN) {
                indexRows = numRows - chunk;
                bTooMany = true;
            }
            for (j = 0; j < indexRows; j++) {
                for (k = 0; k < uR; k++) {
                    permuteMatrix(chunk + j, k) = myCombs(i, indexMatrix(j, k));
                }
            }
            if (bTooMany) {break;}
            chunk += indexRows;
        }
    }
    return permuteMatrix;
}

template <typename typeRcpp, typename typeStd>
typeRcpp MultisetCombination(int n, int r, std::vector<typeStd> v,
                             std::vector<int> Reps, int numRows, bool xtraCol) {
    
    int i, j, numIter, numCols;
    int count = 0, zExpSize = 0;
    // sort v and order Reps by the ordering of v.
    for (i = 0; i < (n-1); i++) {
        for (j = (i+1); j < n; j++) {
            if (v[i] > v[j]) {
                std::swap(v[i], v[j]);
                std::swap(Reps[i], Reps[j]);
            }
        }
    }
    
    std::vector<int> z, zExpand, zIndex, zGroup(r);
    int r1 = r - 1, r2 = r - 2, k = 0;
    for (i = 0; i < n; i++) {zExpSize += Reps[i];}
    
    zIndex.reserve(n);
    zExpand.reserve(zExpSize);
    for (i = 0; i < n; i++) {
        zIndex.push_back(k);
        for (j = 0; j < Reps[i]; j++, k++) {
            zExpand.push_back(i);
        }
    }
    
    z.reserve(r);
    for (i = 0; i < r; i++) {
        z.push_back(zExpand[i]);
    }
    
    if (xtraCol) {numCols  = r + 1;} else {numCols = r;}
    typeRcpp combinationMatrix(numRows, numCols);
    
    while (count < numRows) {
        numIter = n - z[r1];
        if (numIter + count > numRows) {numIter = numRows - count;}
        for (i = 0; i < numIter; i++) {
            for (k = 0; k < r; k++) {combinationMatrix(count, k) = v[zExpand[zIndex[z[k]]]];}
            count++;
            z[r1]++;
        }

        for (i = r2; i >= 0; i--) {
            if (zExpand[zIndex[z[i]]] != zExpand[zExpSize - r + i]) {
                z[i]++;
                zGroup[i] = zIndex[z[i]];
                for (k = (i+1); k < r; k++) {
                    zGroup[k] = zGroup[k-1] + 1;
                    z[k] = zExpand[zGroup[k]];
                }
                break;
            }
        }
    }
    
    return combinationMatrix;
}

template <typename typeRcpp, typename typeStd>
typeRcpp MultisetPermutation(int n, int r, std::vector<typeStd> v,
                             std::vector<int> Reps, int numRows,
                             bool xtraCol, bool bUserCap) {

    int i, j = 0, combRows = 0, count = 0, numCols = 0;
    // sort v and order Reps by the ordering of v.
    for (i = 0; i < (n-1); i++) {
        for (j = (i+1); j < n; j++) {
            if (v[i] > v[j]) {
                std::swap(v[i], v[j]);
                std::swap(Reps[i], Reps[j]);
            }
        }
    }
    
    for (i = 0; i < n; i++) {numCols += Reps[i];}
    std::vector<typeStd> permVec;
    std::vector<int> eachRowCount;
    typename std::vector<typeStd>::iterator myIter, permEnd;
    bool retAllPerms = true;
    typeRcpp myCombs;
    
    if (r < numCols) {
        numCols = r;
        retAllPerms = false;
        combRows = (int)MultisetCombRowNum(n,r,Reps);
        myCombs =  MultisetCombination<typeRcpp>(n,r,v,Reps,combRows,false);
        std::vector<int> stdRowVec(r);
        typeStd temp;
        bool keepGoing;
        double testRows = 0, computedRows;
        eachRowCount.reserve(combRows);
        int rowCount, k;
        
        for (i = 0; i < combRows; i++) {
            j = k = 0;
            while (j < r) {
                stdRowVec[j] = k;
                temp = myCombs(i, j);
                keepGoing = true;
                while (keepGoing) {
                    stdRowVec[j] = k;
                    j++;
                    if (j >= r) {
                        keepGoing = false;
                    } else {
                        if (myCombs(i, j) != temp) {
                            keepGoing = false;
                        }
                    }
                }
                k++;
            }
            
            computedRows = NumPermsWithRep(stdRowVec);
            if (computedRows > INT_MAX) {stop("The number of rows cannot exceed 2^31 - 1.");}
            rowCount = (int)computedRows;
            testRows += rowCount;
            eachRowCount.push_back(rowCount);
        }
        
        if (bUserCap) {
            if (numRows > testRows) {numRows = testRows;}
        } else {
            if (testRows > INT_MAX) {stop("The number of rows cannot exceed 2^31 - 1.");}
            numRows = testRows;
        }
    }
    
    if (xtraCol) {numCols++;}
    typeRcpp permuteMatrix(numRows, numCols);
    
    if (retAllPerms) {
        permVec.reserve(numCols);
        for (i = 0; i < n; i++) {
            for (j = 0; j < Reps[i]; j++) {
                permVec.push_back(v[i]);
            }
        }

        for (i = 0; i < numRows; i++) {
            permEnd = permVec.end();
            for (myIter = permVec.begin(), j = 0; myIter < permEnd; myIter++, j++) {
                permuteMatrix(i, j) = *myIter;
            }
            std::next_permutation(permVec.begin(), permEnd);
        }
    } else {
        
        int k, segment = 0;
        bool bTooMany = false;
        permVec.resize(r);
        
        for (i = 0; i < combRows; i++) {
            // populate permVec with next combination
            for (j = 0; j < r; j++) {permVec[j] = myCombs(i, j);}
            
            if (eachRowCount[i] + count > numRows) {
                segment = numRows - count;
                bTooMany = true;
            } else {
                segment = eachRowCount[i];
            }
            
            // // populate permuteMatrix with all permutations
            // // of a particular combination
            for (j = 0; j < segment; j++, count++) {
                permEnd = permVec.end();
                for (myIter = permVec.begin(), k = 0; myIter < permEnd; myIter++, k++) {
                    permuteMatrix(count, k) = *myIter;
                }
                std::next_permutation(permVec.begin(), permEnd);
            }
            if (bTooMany) {break;}
        }
    }
    
    return permuteMatrix;
}

// [[Rcpp::export]]
SEXP CombinatoricsRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, 
                       SEXP f1, SEXP f2, SEXP lim, SEXP numRow,
                       SEXP RIsComb, SEXP RIsFactor,
                       SEXP RKeepRes, SEXP RFreqs) {
    
    int n, m = 0, i, j, m1, m2;
    int lenFreqs = 0, nRows = 0;
    double computedRows, seqEnd, RUserCap;
    bool IsRepetition, IsInteger;
    bool keepRes, IsComb, IsFactor;
    bool SpecialCase, IsCharacter;
    bool IsConstrained, IsMultiset;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqsExpanded;
    std::vector<std::string > vStr;
    
    switch(TYPEOF(Rv)) {
        case INTSXP: {
            IsCharacter = false;
            IsInteger = true;
            break;
        }
        case REALSXP: {
            IsInteger = IsCharacter = false;
            break;
        }
        case STRSXP: {
            IsCharacter = true;
            IsInteger = false;
            break;
        }
        default: {
            stop("Only integers, numerical, character, and factor classes are supported for v");   
        }
    }
    
    IsComb = as<bool>(RIsComb);
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        myReps.push_back(1);
    } else {
        IsMultiset = true;
        switch(TYPEOF(RFreqs)) {
            case REALSXP: {
                myReps = as<std::vector<int> >(RFreqs);
                break;
            }
            case INTSXP: {
                myReps = as<std::vector<int> >(RFreqs);
                break;
            }
            default: {
                stop("freqs must be of type numeric or integer");
            }
        }
        lenFreqs = myReps.size();
        for (i = 0; i < lenFreqs; i++) {
            if (myReps[i] < 1) {stop("each element in freqs must be a positive number");}
            for (j = 0; j < myReps[i]; j++) {freqsExpanded.push_back(i);}
        }
    }
    
    if (Rf_isNull(Rm)) {
        if (IsMultiset) {
            m = freqsExpanded.size();
        } else {
            stop("m and freqs cannot both be NULL");
        }
    } else {
        if (Rf_length(Rm) > 1) {stop("length of m must be 1");}
        switch(TYPEOF(Rm)) {
            case REALSXP: {
                m = as<int>(Rm);
                break;
            }
            case INTSXP: {
                m = as<int>(Rm);
                break;
            }
            default: {
                stop("m must be of type numeric or integer");
            }
        }
    }
    
    if (m < 1) {stop("m must be positive");}
    std::vector<double> rowVec(m);
    if (!Rf_isLogical(Rrepetition)) {stop("repetitions must be a logical value");}
    IsRepetition = as<bool>(Rrepetition);
    IsFactor = as<bool>(RIsFactor);
    if (!Rf_isLogical(RKeepRes)) {stop("keepResults must be a logical value");}
    keepRes = as<bool>(RKeepRes);

    if (IsCharacter) {
        vStr = as<std::vector<std::string > >(Rv);
        n = vStr.size();
    } else {
        if (Rf_length(Rv) == 1) {
            seqEnd = as<double>(Rv);
            if (NumericVector::is_na(seqEnd)) {seqEnd = 1;}
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            IntegerVector vTemp = seq(m1, m2);
            IsInteger = true;
            vNum = as<std::vector<double> >(vTemp);
        } else {
            vNum = as<std::vector<double> >(Rv);
        }

        for (i = (vNum.size() - 1); i >= 0; i--) {
            if (NumericVector::is_na(vNum[i])) {
                vNum.erase(vNum.begin() + i);
            }
        }
        n = vNum.size();
    }
     
    if (IsInteger) {vInt.assign(vNum.begin(), vNum.end());}
    if (IsFactor) {IsCharacter = IsInteger = false;}
    
    if (Rf_isNull(lim)) {
        IsConstrained = false;
    } else {
        if (Rf_isNull(f1)) {
            IsConstrained = false;
        } else {
            if (!Rf_isString(f1)) {stop("constraintFun must be passed as a character");}
            if (Rf_isNull(f2)) {
                IsConstrained = false;
            } else if (IsFactor) {
                IsConstrained = false;
            } else {
                IsConstrained = true;
                if (!Rf_isString(f2)) {stop("comparisonFun must be passed as a character");}
            }
        }
    }

    if (IsMultiset) {
        if (n != lenFreqs) {stop("the length of freqs must equal the length of v");}
        if (IsComb || IsConstrained) {
            if (m > (int)freqsExpanded.size()) {m = freqsExpanded.size();}
            computedRows = MultisetCombRowNum(n, m, myReps);
        } else {
            if (Rf_isNull(Rm)) {
                computedRows = NumPermsWithRep(freqsExpanded);
            } else if (m == (int)freqsExpanded.size()) {
                computedRows = NumPermsWithRep(freqsExpanded);
            } else {
                if (m > (int)freqsExpanded.size()) {m = freqsExpanded.size();}
                computedRows = MultisetCombRowNum(n, m, myReps);
            }
        }
    } else {
        if (IsRepetition) {
            if (IsComb) {
                computedRows = GetRowNum(n, m);
            } else {
                computedRows = std::pow((double)n, (double)m);
            }
        } else {
            if (m > n) {stop("m must be less than or equal to the length of v");}
            if (IsComb) {
                computedRows = nChooseK(n, m);
            } else {
                computedRows = NumPermsNoRep(n, m);
            }
        }
    }
    
    bool bUserRow = false;
    RUserCap = 0;
    if (!Rf_isNull(numRow)) {
        bUserRow = true;
        switch(TYPEOF(numRow)) {
            case REALSXP: {
                RUserCap = as<double>(numRow);
                break;
            }
            case INTSXP: {
                RUserCap = as<double>(numRow);
                break;
            }
            default: {
                stop("rowCap must be of type numeric or integer");
            }
        }
    }
    
    if (RUserCap == 0) {
        if (computedRows > INT_MAX) {stop("The number of rows cannot exceed 2^31 - 1.");}
        RUserCap = nRows = computedRows;
    } else if (RUserCap < 0) {
        stop("The number of rows must be positive");
    } else if (RUserCap > INT_MAX) {
        stop("The number of rows cannot exceed 2^31 - 1.");
    } else {
        nRows = RUserCap;
        if (nRows > computedRows) {
            if (!IsMultiset || IsComb || !IsConstrained) {nRows = computedRows;}
        }
    }

    if (IsCharacter) {
        if (IsMultiset) {
            if (IsComb) {
                return MultisetCombination<CharacterMatrix>(n, m, vStr, myReps, nRows, false);
            } else {
                return MultisetPermutation<CharacterMatrix>(n, m, vStr, myReps, nRows, false, bUserRow);
            }
        } else {
            if (IsComb){
                return ComboGeneral<CharacterMatrix>(n, m, vStr, IsRepetition, nRows, false);
            } else {
                return PermuteGeneral<CharacterMatrix>(n, m, vStr, IsRepetition, nRows, false);
            }
        }
    } else {
        if (IsConstrained) {
            double myLim;

            switch(TYPEOF(lim)) {
                case REALSXP: {
                    myLim = as<double >(lim);
                    break;
                }
                case INTSXP: {
                    myLim = as<double>(lim);
                    break;
                }
                default: {
                    stop("limitConstraints must be of type numeric or integer");
                }
            }

            if (NumericVector::is_na(myLim)) {stop("limitConstraints cannot be NA");}

            std::string mainFun1 = as<std::string >(f1);
            if (mainFun1 != "prod" && mainFun1 != "sum" && mainFun1 != "mean"
                    && mainFun1 != "max" && mainFun1 != "min") {
                stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
            }

            std::string compFun = as<std::string >(f2);
            if (compFun != "<" && compFun != "<=" && compFun != ">"
                    && compFun != ">=" && compFun != "=="
                    && compFun != "=<" && compFun != "=>") {
                stop("comparisonFun must be one of the following: >, >=, <, <=, or ==");
            }

            if (compFun == "=<") {compFun = "<=";}
            if (compFun == "=>") {compFun = ">=";}

            SpecialCase = false;
            if (mainFun1 == "prod") {
                for (i = 0; i < n; i++) {
                    if (vNum[i] < 0) {
                        SpecialCase = true;
                        break;
                    }
                }
            }

            XPtr<funcPtr> xpFun1 = putFunPtrInXPtr(mainFun1);
            funcPtr myFun1 = *xpFun1;

            if (SpecialCase) {
                if (computedRows > INT_MAX) {stop("The number of rows cannot exceed 2^31 - 1.");}
                nRows = computedRows;
                NumericMatrix matRes;
                double testVal;
                bool Success;
                std::vector<int> indexMatch;
                std::vector<double> matchVals;
                indexMatch.reserve(nRows);
                matchVals.reserve(nRows);
                
                if (IsComb) {
                    if (IsMultiset) {
                        matRes = MultisetCombination<NumericMatrix>(n, m, vNum, myReps, nRows, keepRes);
                    } else {
                        matRes = ComboGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes);
                    }
                } else {
                    if (IsMultiset) {
                        matRes = MultisetPermutation<NumericMatrix>(n, m, vNum, myReps, nRows, keepRes, bUserRow);
                        nRows = matRes.nrow();
                    } else {
                        matRes = PermuteGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes);
                    }
                }
                
                XPtr<compPtr> xpComp = putCompPtrInXPtr(compFun);
                compPtr myComp = *xpComp;

                for (i = 0; i < nRows; i++) {
                    for (j = 0; j < m; j++) {rowVec[j] = matRes(i, j);}
                    testVal = myFun1(rowVec);
                    Success = myComp(testVal, myLim);
                    if (Success) {
                        indexMatch.push_back(i);
                        matchVals.push_back(testVal);
                    }
                }
                
                int numCols = m;
                if (keepRes) {numCols++;}
                computedRows = indexMatch.size();
                if (computedRows > RUserCap) {nRows = RUserCap;} else {nRows = computedRows;}
                NumericMatrix returnMatrix(nRows, numCols);
                
                for (i = 0; i < nRows; i++) {
                    for (j = 0; j < m; j++) {
                        returnMatrix(i,j) = matRes(indexMatch[i],j);
                    }
                    if (keepRes) {returnMatrix(i,m) = matchVals[i];}
                }

                return returnMatrix;
            }

            if (keepRes) {
                NumericMatrix matRes = CombinatoricsConstraints<NumericMatrix>(n, m, vNum, IsRepetition,
                                                                        mainFun1, compFun, myLim, nRows,
                                                                        IsComb, true, myReps, IsMultiset);
                nRows = matRes.nrow();

                for (i = 0; i < nRows; i++) {
                    for (j = 0; j < m; j++) {rowVec[j] = matRes(i, j);}
                    matRes(i, m) = myFun1(rowVec);
                }

                return matRes;
            } else {
                if (IsInteger) {
                    return CombinatoricsConstraints<IntegerMatrix>(n, m, vNum, IsRepetition,
                                                            mainFun1, compFun, myLim, nRows,
                                                            IsComb, false, myReps, IsMultiset);
                } else {
                    return CombinatoricsConstraints<NumericMatrix>(n, m, vNum, IsRepetition,
                                                            mainFun1, compFun, myLim, nRows,
                                                            IsComb, false, myReps, IsMultiset);
                }
            }
        } else {
            if (Rf_isNull(f1)) {keepRes = false;}
            if (keepRes) {
                NumericMatrix matRes;

                std::string mainFun2 = as<std::string >(f1);
                if (mainFun2 != "prod" && mainFun2 != "sum" && mainFun2 != "mean"
                        && mainFun2 != "max" && mainFun2 != "min") {
                    stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
                }

                XPtr<funcPtr> xpFun2 = putFunPtrInXPtr(mainFun2);
                funcPtr myFun2 = *xpFun2;

                if (IsComb) {
                    if (IsMultiset) {
                        matRes = MultisetCombination<NumericMatrix>(n, m, vNum, myReps, nRows, true);
                    } else {
                        matRes = ComboGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, true);
                    }
                } else {
                    if (IsMultiset) {
                        matRes = MultisetPermutation<NumericMatrix>(n, m, vNum, myReps, nRows, true, bUserRow);
                        nRows = matRes.nrow();
                    } else {
                        matRes = PermuteGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, true);
                    }
                }

                for (i = 0; i < nRows; i++) {
                    for (j = 0; j < m; j++) {rowVec[j] = matRes(i, j);}
                    matRes(i, m) = myFun2(rowVec);
                }

                return matRes;
            } else {
                if (IsFactor) {
                    IntegerMatrix factorMat;
                    IntegerVector testFactor = as<IntegerVector>(Rv);
                    CharacterVector myClass = testFactor.attr("class");
                    CharacterVector myLevels = testFactor.attr("levels");

                    if (IsMultiset) {
                        if (IsComb) {
                            factorMat = MultisetCombination<IntegerMatrix>(n, m, vInt, myReps, nRows, false);
                        } else {
                            factorMat = MultisetPermutation<IntegerMatrix>(n, m, vInt, myReps, nRows, false, bUserRow);
                        }
                    } else {
                        if (IsComb) {
                            factorMat = ComboGeneral<IntegerMatrix>(n, m, vInt, IsRepetition, nRows, false);
                        } else {
                            factorMat = PermuteGeneral<IntegerMatrix>(n, m, vInt, IsRepetition, nRows, false);
                        }
                    }
                    
                    factorMat.attr("class") = myClass;
                    factorMat.attr("levels") = myLevels;
                    
                    return factorMat;
                } else {
                    if (IsComb) {
                        if (IsInteger) {
                            if (IsMultiset) {
                                return MultisetCombination<IntegerMatrix>(n, m, vInt, myReps, nRows, false);
                            } else {
                                return ComboGeneral<IntegerMatrix>(n, m, vInt, IsRepetition, nRows, false);
                            }
                        } else {
                            if (IsMultiset) {
                                return MultisetCombination<NumericMatrix>(n, m, vNum, myReps, nRows, false);
                            } else {
                                return ComboGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, false);
                            }
                        }
                    } else {
                        if (IsInteger) {
                            if (IsMultiset) {
                                return MultisetPermutation<IntegerMatrix>(n, m, vInt, myReps, nRows, false, bUserRow);
                            } else {
                                return PermuteGeneral<IntegerMatrix>(n, m, vInt, IsRepetition, nRows, false);
                            }
                        } else {
                            if (IsMultiset) {
                                return MultisetPermutation<NumericMatrix>(n, m, vNum, myReps, nRows, false, bUserRow);
                            } else {
                                return PermuteGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, false);
                            }
                        }
                    }
                }
            }
        }
    }
}
