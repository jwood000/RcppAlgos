#include <Rcpp.h>
#include <algorithm>
#include <string>
using namespace Rcpp;

List rleCpp(std::vector<double> x) {
    std::vector<unsigned long int> lengths, numUni;
    std::vector<double> values;
    std::vector<double>::iterator it, xBeg, xEnd;
    xBeg = x.begin() + 1; xEnd = x.end();
    double prev = x[0];
    unsigned long int n = x.size(), i = 0;
    lengths.reserve(n);
    values.reserve(n);
    values.push_back(prev);
    lengths.push_back(1);
    for(it = xBeg; it < xEnd; it++) {
        if (prev == *it) {
            lengths[i]++;
        } else {
            values.push_back(*it);
            lengths.push_back(1);
            i++;
            prev = *it;
        }
    }
    
    numUni.push_back(i);
    
    return List::create(
        _["lengths"] = lengths,
        _["values"] = values,
        _["uniques"] = numUni
    );
}

int NumPermsWithRep(std::vector<double> v) {
    List myRle = rleCpp(v);
    unsigned long int n = v.size(), myMax;
    std::vector<unsigned long int> myLens = myRle[0], myUnis = myRle[2];
    std::sort(myLens.begin(), myLens.end(),
              std::greater<unsigned long int>());
    
    myMax = myLens[0];
    unsigned long int i, j, numUni = myUnis[0];
    int result = 1;
    
    for (i = n; i > myMax; i--) {result *= i;}

    if (numUni > 0) {
        for (i = 1; i <= numUni; i++) {
            // No need to divide by 1.
            // Start j at 2 instead.
            for (j = 2; j <= myLens[i]; j++) {
                result /= j;
            }
        }
    }

    return result;
}

double NumPermsNoRep(int n, int k) {
    double dblN = (double)n, result = 1;
    double i, m = dblN - (double)k;
    for (i = n; i > m; i--) {result *= i;}
    return result;
}

template <typename TypeRcpp>
TypeRcpp SubMat(TypeRcpp m, int n) {
    int k = m.ncol();
    TypeRcpp subMatrix(n,k);
    for (int i = 0; i < n; i++) {subMatrix(i,_) = m(i,_);}
    return subMatrix;
}

double nChooseK(double n, double k) {
    // returns the number of k-combinations from a set
    // of n elements. Mathematically speaking, 
    //  we have: n!/(k!*(n-k)!)
    double nCk;
    double temp = 1;
    for(int i = 1; i <= k; i++) {temp *= (n - k + i)/i;}
    nCk = round(temp);
    return nCk;
}

double GetRowNum(int n, int r) {
    // for combinations where repetition is allowed, this
    // function returns the number of combinations for
    // a given n and r. The resulting vector, "triangleVec"
    // resembles triangle numbers. In fact, this vector
    // is obtained in a very similar method as generating
    // triangle numbers, albeit in a repeating fashion.
    int i, k;
    std::vector<double> triangleVec(n);
    std::vector<double> temp;
    for (i = 0; i < n; i++) {triangleVec[i] = i+1;}
    
    for (i = 1; i < r; i++) {
        temp.clear();
        temp.reserve(n);
        for (k = 1; k <= n; k++) {
            temp.push_back(std::accumulate(triangleVec.begin(), triangleVec.begin() + k, 0.0));
        }
        triangleVec = temp;
    }
    
    return triangleVec[n-1];
}

// Below, we define five functions that will be utilized
// as constraint functions. We also define five comparison
// operations (<, <=, >, >=, ==). The framework is based on
// the information posted by Dirk Eddelbuettel from
// this Rcpp Gallery (http://gallery.rcpp.org/articles/passing-cpp-function-pointers/)

double prodCpp(std::vector<double>& v, int& r) {
    std::vector<double>::iterator i;
    std::vector<double>::iterator myEnd = v.begin()+r;
    double myProduct = 1.0;
    for (i = v.begin(); i < myEnd; i++) {myProduct *= *i;}
    return(myProduct);
}

double sumCpp(std::vector<double>& v, int& r) {
    std::vector<double>::iterator i;
    std::vector<double>::iterator myEnd = v.begin()+r;
    double mySum = 0.0;
    for (i = v.begin(); i < myEnd; i++) {mySum += *i;}
    return(mySum);
}

double meanCpp(std::vector<double>& v, int& r){
    double mySum = sumCpp(v, r);
    return mySum/r;
}

double maxCpp(std::vector<double>& v, int& r) {
    std::vector<double>::iterator myEnd = v.begin()+r;
    std::vector<double>::iterator y = std::max_element(v.begin(), myEnd);
    return v[std::distance(v.begin(), y)];
}

double minCpp(std::vector<double>& v, int& r) {
    std::vector<double>::iterator myEnd = v.begin()+r;
    std::vector<double>::iterator y = std::min_element(v.begin(), myEnd);
    return v[std::distance(v.begin(), y)];
}

// Comparison functions
bool lessCpp(double& x, double& y) {return x < y;}
bool greaterCpp(double& x, double& y) {return x > y;}
bool lessEqualCpp(double& x, double& y) {return x <= y;}
bool greaterEqualCpp(double& x, double& y) {return x >= y;}
bool equalCpp(double& x, double& y) {return x == y;}

typedef double (*funcPtr)(std::vector<double>& x, int& y);
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

IntegerMatrix MakeIndexHeaps(unsigned long int indRows, unsigned long int r) {
    unsigned long int j, i = 0, count = 0;
    IntegerMatrix indexMatrix(indRows, r);
    
    std::vector<unsigned long int> vecInd(r, 0);
    IntegerVector mySeq = seq(0, r-1);
    for (j = 0; j < r; j++) {indexMatrix(count, j) = j;}
    
    while (i < r) {
        if (vecInd[i] < i) {
            if (i % 2 == 0)
                std::swap(mySeq[0], mySeq[i]);
            else
                std::swap(mySeq[vecInd[i]], mySeq[i]);
            count++;
            for (j = 0; j < r; j++) {indexMatrix(count,j) = mySeq[j];}
            vecInd[i]++;
            i = 0;
        } else {
            vecInd[i] = 0;
            i++;
        }
    }
    
    return indexMatrix;
}

// This function applys a constraint function to a vector v with respect
// to a constraint value "lim". The main idea is that combinations are added
// successively, until a particular combination exceeds the given constraint
// value for a given constraint function. After this point, we can safely skip
// several combinations knowing that they will exceed the given constraint value.
template <typename TypeRcpp>
TypeRcpp CombinatoricsConstraints(int n, int r, std::vector<double> v,
                                       bool repetition, std::string myFun, 
                                       std::string myComparison, double lim,
                                       int rowNum, bool isComb, bool xtraCol) {
    // where myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max"
    // and myComparison is a comparison operator: "<", "<=", ">", or ">="
    
    double extremeV, testVal;
    int count = 0, s = v.size();
    unsigned int ur = r, numCols;
    
    // constraintFun is a pointer to one of the functions defined above
    XPtr<funcPtr> xpFun = putFunPtrInXPtr(myFun);
    funcPtr constraintFun = *xpFun;
    
    // comparisonFun is a pointer to one of the comparisons defined above
    XPtr<compPtr> xpCompOne = putCompPtrInXPtr(myComparison);
    compPtr comparisonFunOne = *xpCompOne;
    
    XPtr<compPtr> xpCompTwo = xpCompOne;
    compPtr comparisonFunTwo;

    if (myComparison == ">" || myComparison == ">=") {
        std::sort(v.begin(), v.end(), std::greater<double>());
        extremeV = minCpp(v, s);
        comparisonFunTwo = *xpCompOne;
    } else {
        std::sort(v.begin(), v.end());
        extremeV = maxCpp(v, s);
        if (myComparison == "==") {
            XPtr<compPtr> xpCompThree = putCompPtrInXPtr("<=");
            comparisonFunTwo = *xpCompThree;
        } else {
            comparisonFunTwo = *xpCompOne;
        }
    }
    
    if (xtraCol) {numCols = r+1;} else {numCols = r;}
    TypeRcpp combinatoricsMatrix(rowNum, numCols);
    std::vector<double> z, zPart;
    bool t_1, t_2, t = true, keepGoing = true;
    std::vector<double>::iterator it;
    std::vector<double> temp;
    int posMaxZ, tSize = r - 1, newSize;
    int numPerms, rMinus, k, i, j, zSize;
    
    if (repetition) {
        v.erase(std::unique(v.begin(), v.end()), v.end());
        s = v.size();
        z.assign(r+s-1, v[0]);
        std::copy(v.begin()+1, v.end(), z.begin()+r);
        while (keepGoing) {
            t_2 = z.size() >= ur;
            testVal = constraintFun(z, r);
            t = comparisonFunTwo(testVal, lim);
            while (t && t_2 && keepGoing) {
                testVal = constraintFun(z, r);
                t_1 = comparisonFunOne(testVal, lim);
                if (t_1) {
                    if (isComb) {
                        for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = z[k];}
                        count++;
                    } else {
                        zPart.clear(); zPart.reserve(r);
                        for (k=0; k < r; k++) {zPart.push_back(z[k]);}
                        numPerms = NumPermsWithRep(zPart);
                        for (i=0; i < numPerms; i++) {
                            for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = zPart[k];}
                            std::next_permutation(zPart.begin(), zPart.end());
                            count++;
                            if (count >= rowNum) {break;}
                        }
                    }
                    keepGoing = (count < rowNum);
                }
                zSize = z.size();
                t_2 = zSize >= r;
                if (t_2) {
                    if (zSize==r) t_2 = false;
                    z.erase(z.begin() + r - 1);
                    testVal = constraintFun(z, r);
                    t = comparisonFunTwo(testVal, lim);
                }
            }
            
            if (keepGoing) {
                if (t) {
                    posMaxZ = -1;
                    zSize = z.size();
                    for (k = zSize-1; k >= 0; k--) {
                        if (z[k] != extremeV) {
                            posMaxZ = k;
                            break;
                        }
                    }
        
                    if (posMaxZ > -1) {
                        for (k = 0; k < s; k++) {
                            if (v[k] == z[posMaxZ]) {
                                z[posMaxZ] = v[k+1];
                                break;
                            }
                        }

                        temp.clear();
                        temp.reserve(posMaxZ);
                        for (it = z.begin(); it <= z.begin() + posMaxZ; it++) {temp.push_back(*it);}
                        
                        if (z[posMaxZ] != extremeV) {
                            newSize = tSize + s - k - 1;
                            z.assign(newSize, z[posMaxZ]);
                            std::copy(temp.begin(), temp.end(), z.begin());
                            std::copy(v.begin() + (k+1), v.end(), z.begin() + tSize);
                        } else {
                            z.assign(r, z[posMaxZ]);
                            std::copy(temp.begin(), temp.end(), z.begin());
                        }
                    } else {
                        keepGoing = false;
                    }
                } else {
                    rMinus = r - 2;
                    while (rMinus >= 0 && !t) {
                        if (z[rMinus] != extremeV) {
                            for (k = 0; k < s; k++) {
                                if (v[k] == z[rMinus]) {
                                    z[rMinus] = v[k+1];
                                    break;
                                }
                            }
                            temp.clear();
                            for (it = z.begin(); it < z.begin() + rMinus; it++) {temp.push_back(*it);}
                            newSize = tSize + s - k - 1;
                            z.assign(newSize, z[rMinus]);
                            std::copy(temp.begin(), temp.end(), z.begin());
                            std::copy(v.begin() + (k+1), v.end(), z.begin() + tSize);
                            if (z.size() >= ur) {
                                testVal = constraintFun(z, r);
                                t = comparisonFunTwo(testVal, lim);
                            } else {t = false;}
                        }
                        --rMinus;
                    }
                    if (!t) {keepGoing = false;}
                }
            }
        }
    } else {
        z = v;
        int indexRows;
        IntegerMatrix indexMatrix;
        if (!isComb) {
            indexRows = (int)NumPermsNoRep(r, r-1);
            indexMatrix = MakeIndexHeaps(indexRows, r);
        }

        while (keepGoing) {
            t_2 = z.size() >= ur;
            testVal = constraintFun(z, r);
            t = comparisonFunTwo(testVal, lim);
            while (t && t_2 && keepGoing) {
                testVal = constraintFun(z, r);
                t_1 = comparisonFunOne(testVal, lim);
                if (t_1) {
                    if (isComb) {
                        for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = z[k];}
                        count++;
                    } else {
                        for (j = 0; j < indexRows; j++) {
                            for (k = 0; k < r; k++) {
                                combinatoricsMatrix(count, k) = z[indexMatrix(j, k)];
                            }
                            count++;
                            if (count >= rowNum) {break;}
                        }
                    }
                    keepGoing = (count < rowNum);
                }
                zSize = z.size();
                t_2 = zSize >= r;
                if (t_2) {
                    if (zSize==r) t_2 = false;
                    z.erase(z.begin() + r - 1);
                    testVal = constraintFun(z, r);
                    t = comparisonFunTwo(testVal, lim);
                }
            }
            
            if (keepGoing) {
                if (t) {
                    zSize = z.size();
                    posMaxZ = -1;
                    for (k = zSize-1; k >= 0; k--) {
                        if (z[k] != v[(s-r+k)]) {
                            posMaxZ = k;
                            break;
                        }
                    }
                    
                    if (posMaxZ > -1) {
                        for (k = 0; k < s; k++) {
                            if (v[k] == z[posMaxZ]) {
                                z[posMaxZ] = v[k+1];
                                break;
                            }
                        }
                        
                        temp.clear();
                        temp.reserve(posMaxZ + s - k - 2);
                        for (it = z.begin(); it <= z.begin() + posMaxZ; it++) {temp.push_back(*it);}
                        for (it = v.begin() + (k+2); it < v.end(); it++) {temp.push_back(*it);}
                        z = temp;
                    } else {
                        keepGoing = false;
                    }
                } else {
                    rMinus = r - 2;
                    while (rMinus >= 0 && !t) {
                        for (k = 0; k < s; k++) {
                            if (v[k] == z[rMinus]) {
                                z[rMinus] = v[k+1];
                                break;
                            }
                        }
                        
                        if (z[rMinus] != extremeV) {
                            temp.clear();
                            temp.reserve(rMinus + s - k - 2);
                            for (it = z.begin(); it <= z.begin() + rMinus; it++) {temp.push_back(*it);}
                            for (it = v.begin() + (k+2); it < v.end(); it++) {temp.push_back(*it);}
                            z = temp;
                            if (z.size() >= ur) {
                                testVal = constraintFun(z, r);
                                t = comparisonFunTwo(testVal, lim);
                            } else {t = false;}
                        }
                        --rMinus;
                    }
                    if (!t) {keepGoing = false;}
                }
            }
        }
    }
    return(SubMat(combinatoricsMatrix, count));
}

IntegerMatrix ComboFactor(int n, int r, IntegerVector v, bool repetition, int rowNum) {
    v.sort();
    IntegerVector z;
    int count = 0, tSize = r - 1, s = v.size();
    int k, i, numIter, posMaxZ, newSize;
    int maxV = v[s-1];
    IntegerVector temp;
    bool keepGoing = true;
    IntegerMatrix combinationMatrix(rowNum, r);
    combinationMatrix.attr("levels") = v.attr("levels");
    combinationMatrix.attr("class") = v.attr("class");
    
    if (repetition) {
        v = unique(v);
        v.sort();
        s = v.size();
        newSize = r + s - 1;
        z = rep(v[0], newSize);
        for (i = r; i < newSize; i++) {z[i] = v[i-r+1];}
        
        while (keepGoing) {
            numIter = z.size() - r;
            for (i = 0; i < numIter; i++) {
                for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
                count++;
                z.erase(z.begin() + tSize);
            }
            
            for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
            count++;
            
            posMaxZ = -1;
            for (k = z.size()-1; k >= 0; k--) {
                if (z[k] < maxV) {
                    posMaxZ = k;
                    break;
                }
            }
            
            if (posMaxZ > -1) {
                for (k = 0; k < s; k++) {
                    if (v[k] == z[posMaxZ]) {
                        z[posMaxZ] = v[k+1];
                        break;
                    }
                }
                
                temp.assign(z.begin(), z.begin() + posMaxZ + 1);
                
                if (z[posMaxZ] < maxV) {
                    newSize = tSize + s - k - 1;
                    z = rep(z[posMaxZ], newSize);
                    std::copy(temp.begin(), temp.end(), z.begin());
                    std::copy(v.begin() + (k+1), v.end(), z.begin() + tSize);
                } else {
                    z = rep(z[posMaxZ], r);
                    std::copy(temp.begin(), temp.end(), z.begin());
                }
            } else {
                keepGoing = false;
            }
        }
    } else {
        z = v;
        while (keepGoing) {
            numIter = z.size() - r;
            for (i = 0; i < numIter; i++) {
                for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
                count++;
                z.erase(z.begin() + tSize);
            }
            
            for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
            count++;
            
            posMaxZ = -1;
            for (k = z.size()-1; k >= 0; k--) {
                if (z[k] != v[s-r+k]) {
                    posMaxZ = k;
                    break;
                }
            }
            
            if (posMaxZ > -1) {
                for (k = 0; k < s; k++) {
                    if (v[k] == z[posMaxZ]) {
                        z[posMaxZ] = v[k+1];
                        break;
                    }
                }
                
                temp.assign(z.begin(), z.begin() + posMaxZ + 1);
                for (i = (k+2); i < s; i++) {temp.push_back(v[i]);}
                z = temp;
            } else {
                keepGoing = false;
            }
        }
    }
    return(combinationMatrix);
}

IntegerMatrix PermuteFactor(int n, int r, IntegerVector v, bool repetition, int rowNum) {
    unsigned long int uN = n, uR = r, uRowN = rowNum;
    unsigned long int i = 0, j, k, chunk;
    IntegerMatrix permuteMatrix(uRowN, uR);
    permuteMatrix.attr("class") = v.attr("class");
    permuteMatrix.attr("levels") = v.attr("levels");
    
    if (repetition) {
        unsigned long int groupLen = 1, repLen = 1;
        IntegerVector::iterator m, vBeg, vEnd;
        vBeg = v.begin(); vEnd = v.end();
        for (i = 0; i < uR; i++) {
            groupLen *= uN;
            chunk = 0;
            for (k = 0; k < uRowN; k += groupLen) {
                for (m = vBeg; m < vEnd; m++) {
                    for (j = 0; j < repLen; j++) {
                        permuteMatrix(chunk + j, i) = *m;
                    }
                    chunk += repLen;
                }
            }
            repLen *= uN;
        }
    } else {
        unsigned long int combRows = (unsigned long int)nChooseK(n, r);
        IntegerMatrix myCombs = ComboFactor(n,r,v,false,combRows);
        unsigned long int indexRows = (unsigned long int)NumPermsNoRep(r, r-1);
        IntegerMatrix indexMatrix = MakeIndexHeaps(indexRows, uR);
        
        chunk = 0;
        for (i = 0; i < combRows; i++) {
            for (j = 0; j < indexRows; j++) {
                for (k = 0; k < uR; k++) {
                    permuteMatrix(chunk + j, k) = myCombs(i, indexMatrix(j, k));
                }
            }
            chunk += indexRows;
        }
    }
    return(permuteMatrix);
}

template <typename TypeRcpp, typename stdType>
TypeRcpp ComboGeneral(int n, int r, std::vector<stdType> v,
                      bool repetition, int rowNum, bool xtraCol) {
    std::sort(v.begin(), v.end());
    std::vector<stdType> z;
    int count = 0, tSize = r - 1, s = v.size();
    int k, i, numIter, posMaxZ, newSize, numCols;
    stdType maxV = v[s-1];
    typename std::vector<stdType>::iterator it;
    std::vector<stdType> temp;
    bool keepGoing = true;
    if (xtraCol) {numCols  = r + 1;} else {numCols = r;}
    TypeRcpp combinationMatrix(rowNum, numCols);
    
    if (repetition) {
        v.erase(std::unique(v.begin(), v.end()), v.end());
        s = v.size();
        z.assign(r+s-1, v[0]);
        std::copy(v.begin()+1, v.end(), z.begin()+r);
        
        while (keepGoing) {
            numIter = z.size() - r;
            for (i = 0; i < numIter; i++) {
                for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
                count++;
                z.erase(z.begin() + tSize);
            }
            
            for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
            count++;
            
            posMaxZ = -1;
            for (k = z.size()-1; k >= 0; k--) {
                if (z[k] < maxV) {
                    posMaxZ = k;
                    break;
                }
            }
            
            if (posMaxZ > -1) {
                for (k = 0; k < s; k++) {
                    if (v[k] == z[posMaxZ]) {
                        z[posMaxZ] = v[k+1];
                        break;
                    }
                }
                
                temp.clear();
                temp.reserve(posMaxZ);
                for (it = z.begin(); it <= z.begin() + posMaxZ; it++) {temp.push_back(*it);}
                
                if (z[posMaxZ] < maxV) {
                    newSize = tSize + s - k - 1;
                    z.assign(newSize, z[posMaxZ]);
                    std::copy(temp.begin(), temp.end(), z.begin());
                    std::copy(v.begin() + (k+1), v.end(), z.begin() + tSize);
                } else {
                    z.assign(r, z[posMaxZ]);
                    std::copy(temp.begin(), temp.end(), z.begin());
                }
            } else {
                keepGoing = false;
            }
        }
    } else {
        z = v;
        while (keepGoing) {
            numIter = z.size() - r;
            for (i = 0; i < numIter; i++) {
                for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
                count++;
                z.erase(z.begin() + tSize);
            }
            
            for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
            count++;
            
            posMaxZ = -1;
            for (k = z.size()-1; k >= 0; k--) {
                if (z[k] != v[(s-r+k)]) {
                    posMaxZ = k;
                    break;
                }
            }
            
            if (posMaxZ > -1) {
                for (k = 0; k < s; k++) {
                    if (v[k] == z[posMaxZ]) {
                        z[posMaxZ] = v[k+1];
                        break;
                    }
                }
                
                temp.clear();
                temp.reserve(posMaxZ + s - k - 2);
                for (it = z.begin(); it <= z.begin() + posMaxZ; it++) {temp.push_back(*it);}
                for (it = v.begin() + (k+2); it < v.end(); it++) {temp.push_back(*it);}
                z = temp;
            } else {
                keepGoing = false;
            }
        }
    }
    return(combinationMatrix);
}

template <typename TypeRcpp, typename stdType>
TypeRcpp PermuteGeneral(int n, int r, std::vector<stdType> v,
                        bool repetition, int rowNum, bool xtraCol) {
    unsigned long int uN = n, uR = r, uRowN = rowNum;
    unsigned long int i = 0, j, k, chunk, numCols;
    if (xtraCol) {numCols  = uR + 1;} else {numCols = uR;}
    TypeRcpp permuteMatrix(uRowN, numCols);
    
    if (repetition) {
        unsigned long int groupLen = 1, repLen = 1;
        typename std::vector<stdType>::iterator m, vBeg, vEnd;
        vBeg = v.begin(); vEnd = v.end();
        for (i = 0; i < uR; i++) {
            groupLen *= uN;
            chunk = 0;
            for (k = 0; k < uRowN; k += groupLen) {
                for (m = vBeg; m < vEnd; m++) {
                    for (j = 0; j < repLen; j++) {
                        permuteMatrix(chunk + j, i) = *m;
                    }
                    chunk += repLen;
                }
            }
            repLen *= uN;
        }
    } else {
        unsigned long int combRows = (unsigned long int)nChooseK(n, r);
        TypeRcpp myCombs = ComboGeneral<TypeRcpp>(n,r,v,false,combRows,false);
        unsigned long int indexRows = (unsigned long int)NumPermsNoRep(r, r-1);
        IntegerMatrix indexMatrix = MakeIndexHeaps(indexRows, uR);

        chunk = 0;
        for (i = 0; i < combRows; i++) {
            for (j = 0; j < indexRows; j++) {
                for (k = 0; k < uR; k++) {
                    permuteMatrix(chunk + j, k) = myCombs(i, indexMatrix(j, k));
                }
            }
            chunk += indexRows;
        }
    }
    return(permuteMatrix);
}

int SectionLength2(std::vector<int> v, unsigned long int n) {
    unsigned long int numUni = v.size();
    std::sort(v.begin(), v.end(),
              std::greater<unsigned long int>());
    
    unsigned long int i, j, myMax;
    int result = 1;
    
    myMax = v[0];
    for (i = n; i > myMax; i--) {result *= i;}
    
    if (numUni > 1) {
        for (i = 1; i < numUni; i++) {
            // No need to divide by 1.
            // Start j at 2 instead.
            for (j = 2; j <= v[i]; j++) {
                result /= j;
            }
        }
    }
    
    return result;
}

std::vector<unsigned long int> SectionLength(std::vector<int> v,
                                             std::vector<int> repLen) {
    unsigned long int n = v.size(), mySum;
    unsigned long int i, j, k, numRow = 1;
    unsigned long int n2 = n-2, n1 = n-1, res;
    std::vector<int>::iterator it;
    std::vector<int> temp(n);
    for (it = v.begin(); it < v.end(); it++) {numRow *= (*it+1);}
    std::vector<unsigned long int> myLengths(numRow);
    
    for (i = 0; i < numRow; i++) {
        temp[0] = i / repLen[0];
        mySum = temp[0];
        temp[n1] = (i % repLen[n2]);
        mySum += temp[n1];
        for (j = 1; j < n1; j++) {
            temp[j] = (i % repLen[j-1]) / repLen[j];
            mySum += temp[j];
        }
    
        std::sort(temp.begin(), temp.end(),
                  std::greater<unsigned long int>());
        res = 1;
        for (j = mySum; j > temp[0]; j--) {res *= j;}
        
        for (j = 1; j < n; j++) {
            for (k = 2; k <= temp[j]; k++) {
                res /= k;
            }
        }
    
        myLengths[i] = res;
    }

    return myLengths;
}

//template <typename TypeRcpp, typename stdType>
// PermuteSpecial(3, c(0,1,2), c(4,2,2), 420)
// [[Rcpp::export]]
IntegerMatrix PermuteSpecial(int n, std::vector<int> v,
                             std::vector<int> Reps, int rowNum) {
    unsigned long int uRowN = rowNum, uN = n, count;
    unsigned long int i, j, k, m = 1, r, numCols, strt, ind;

    numCols = std::accumulate(Reps.begin(), Reps.end(), 0);
    unsigned long int lastCol = numCols - 1, lastSum = 0;
    IntegerMatrix permuteMatrix(uRowN, numCols);

    std::vector<int> specPerm(uN);
    std::vector<unsigned long int> vecLast(uRowN);
    
    int myInt;
    std::vector<int> repLen(n, 1);
    for (myInt = (n-2); myInt >= 0; myInt--) {
        repLen[myInt] = repLen[myInt+1]*(Reps[myInt+1]+1);
    }
    
    std::vector<unsigned long int> groupLen = SectionLength(Reps, repLen);
    
    std::vector<std::vector<int> > prevPerms(2*uRowN/5, std::vector<int>(uN));
    std::vector<std::vector<int> > nextPerms(4*uRowN/5, std::vector<int>(uN));
    std::vector<int>::iterator it, vBeg, vEnd;
    vBeg = v.begin(); vEnd = v.end();
    numCols--;
    nextPerms[0] = Reps;
    
    for (i = 1; i < uN; i++) {
        for (j = 0; j < Reps[i]; j++) {
            lastSum += i;
        }
    }

    for (i = 0; i < uRowN; i++) {vecLast[i] = lastSum;}

    for (i = 0; i < (lastCol-1); i++) {
        for (j = 0; j < m; j++) {prevPerms[j] = nextPerms[j];}
        strt = count = m = r = 0;
        while (count < uRowN) {
            specPerm = prevPerms[r];
            j = ind = 0;
            for (k = 0; k < uN; k++) {ind += specPerm[k]*repLen[k];}
            for (it = vBeg; it < vEnd; it++) {
                if (specPerm[j] > 0) {
                    specPerm[j]--;
                    ind -= repLen[j];
                    count += groupLen[ind];
                    for (k = strt; k < count; k++) {
                        permuteMatrix(k, i) = *it;
                        vecLast[k] -= j;
                    }
                    strt = count;
                    nextPerms[m] = specPerm;
                    m++;
                    specPerm[j]++;
                    ind += repLen[j];
                }
                j++;
            }
            r++;
        }
    }
    
    for (i = (lastCol-1); i < lastCol; i++) {
        strt = count = m = r = 0;
        while (count < uRowN) {
            specPerm = nextPerms[r];
            j = ind = 0;
            for (k = 0; k < uN; k++) {ind += specPerm[k]*repLen[k];}
            for (it = vBeg; it < vEnd; it++) {
                if (specPerm[j] > 0) {
                    specPerm[j]--;
                    ind -= repLen[j];
                    count += groupLen[ind];
                    for (k = strt; k < count; k++) {
                        permuteMatrix(k, i) = *it;
                        vecLast[k] -= j;
                    }
                    strt = count;
                    specPerm[j]++;
                    ind += repLen[j];
                }
                j++;
            }
            r++;
        }
    }

    for (i = 0; i < uRowN; i++) {
        permuteMatrix(i, lastCol) = v[vecLast[i]];
    }
    
    return(permuteMatrix);
}

// [[Rcpp::export]]
IntegerMatrix PermuteSpecial2(int n, std::vector<int> v,
                             std::vector<int> Reps, int rowNum) {
    unsigned long int uRowN = rowNum, uN = n;
    unsigned long int i, j, k, m, r, numCols, strt;
    unsigned long int groupLen = 1, count;
    
    numCols = std::accumulate(Reps.begin(), Reps.end(), 0);
    unsigned long int lastCol = numCols - 1, lastSum = 0;
    IntegerMatrix permuteMatrix(uRowN, numCols);
    
    std::vector<int> specPerm(uN);
    std::vector<unsigned long int> vecLast(uRowN);
    std::vector<std::vector<int> > prevPerms(uRowN, std::vector<int>(uN));
    std::vector<std::vector<int> > nextPerms(uRowN, std::vector<int>(uN));
    numCols--;
    nextPerms[0] = Reps;
    
    for (i = 1; i < uN; i++) {
        for (j = 0; j < Reps[i]; j++) {
            lastSum += i;
        }
    }
    
    for (i = 0; i < uRowN; i++) {vecLast[i] = lastSum;}
    
    for (i = 0; i < (lastCol-1); i++) {
        prevPerms = nextPerms;
        strt = count = m = r = 0;
        while (count < uRowN) {
            specPerm = prevPerms[r];
            for (j = 0; j < uN; j++) {
                if (specPerm[j] > 0) {
                    specPerm[j]--;
                    groupLen = SectionLength2(specPerm, numCols-i);
                    count += groupLen;
                    for (k = strt; k < count; k++) {
                        permuteMatrix(k, i) = v[j];
                        vecLast[k] -= j;
                    }
                    strt = count;
                    nextPerms[m] = specPerm;
                    m++;
                    specPerm[j]++;
                }
            }
            r++;
        }
    }
    
    for (i = (lastCol-1); i < lastCol; i++) {
        strt = count = m = r = 0;
        while (count < uRowN) {
            specPerm = nextPerms[r];
            for (j = 0; j < uN; j++) {
                if (specPerm[j] > 0) {
                    specPerm[j]--;
                    groupLen = SectionLength2(specPerm, numCols-i);
                    count += groupLen;
                    for (k = strt; k < count; k++) {
                        permuteMatrix(k, i) = v[j];
                        vecLast[k] -= j;
                    }
                    strt = count;
                    m++;
                    specPerm[j]++;
                }
            }
            r++;
        }
    }
    
    for (i = 0; i < uRowN; i++) {
        permuteMatrix(i, lastCol) = v[vecLast[i]];
    }
    
    return(permuteMatrix);
}

// [[Rcpp::export]]
SEXP CombinatoricsRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, 
                       SEXP f1, SEXP f2, SEXP lim, SEXP numRow,
                       SEXP RIsComb, SEXP RIsFactor, SEXP RKeepRes) {
    int n, m, i, j, m1, m2, nRows = 0;
    double testRows;
    bool IsRepetition, IsConstrained;
    bool IsCharacter, IsComb, IsFactor;
    bool SpecialCase, keepRes, IsInteger;
    
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
    
    if (m < 1) {stop("m must be positive");}
    std::vector<double> vNum;
    IntegerVector vFactor;
    std::vector<int> vInt;
    std::vector<std::string > vStr;
    std::vector<double> rowVec(m);
    
    switch(TYPEOF(Rv)) {
        case REALSXP: {
            IsCharacter = false;
            IsInteger = false;
            break;
        }
        case INTSXP: {
            IsCharacter = false;
            IsInteger = true;
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
    
    if (!Rf_isLogical(Rrepetition)) {stop("repetitions must be a logical value");}
    IsRepetition = as<bool >(Rrepetition);
    IsComb = as<bool>(RIsComb);
    IsFactor = as<bool>(RIsFactor);
    keepRes = as<bool>(RKeepRes);
    
    if (IsCharacter) {
        vStr = as<std::vector<std::string > >(Rv);
        n = vStr.size();
    } else {
        if (Rf_length(Rv) == 1) {
            j = as<int>(Rv);
            if (j > 1) {m1 = 1; m2 = j;} else {m1 = j; m2 = 1;}
            IntegerVector vTemp = seq(m1, m2);
            IsInteger = true;
            vNum = as<std::vector<double> >(vTemp);
        } else {
            vNum = as<std::vector<double> >(Rv);    
        }
        n = vNum.size();
    }
    
    if (IsInteger) {vInt.assign(vNum.begin(), vNum.end());}
    
    if (IsFactor) {
        IntegerVector testFactor = as<IntegerVector>(Rv);
        vFactor = seq(1, n);
        CharacterVector myClass = testFactor.attr("class");
        CharacterVector myLevels = testFactor.attr("levels");
        vFactor.attr("class") = myClass;
        vFactor.attr("levels") = myLevels;
        IsCharacter = IsInteger = false;
    }
    
    if (IsRepetition) {
        if (IsComb) {
            testRows = GetRowNum(n, m);
        } else {
            testRows = pow((double)n, (double)m);
        }
    } else {
        if (m > n) {stop("m must be less than or equal to the length of v");}
        if (IsComb) {
            testRows = nChooseK(n, m);
        } else {
            testRows = NumPermsNoRep(n, m);
        }
    }

    if (IsCharacter) {
        if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
        nRows = testRows;
        if (IsComb){
            return ComboGeneral<CharacterMatrix>(n, m, vStr, IsRepetition, nRows, false);
        } else {
            return PermuteGeneral<CharacterMatrix>(n, m, vStr, IsRepetition, nRows, false);
        }
    } else {
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
        
        if (IsConstrained) {
            double myLim, RUserCap;

            switch(TYPEOF(lim)) {
                case REALSXP: {
                    myLim = as<double >(lim);
                    break;
                }
                case INTSXP: {
                    myLim = as<double >(lim);
                    break;
                }
                default: {
                    stop("limitConstraints must be of type numeric or integer");
                }
            }

            if (Rf_isNull(numRow)) {
                RUserCap = 0;
            } else {
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
            
            if (RUserCap == 0) {
                if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
                RUserCap = nRows = testRows;
            } else if (RUserCap < 0) {
                stop("The number of rows must be positive");
            } else if (RUserCap > 2147483647) {
                stop("The number of rows cannot exceed 2^31 - 1.");
            } else {
                if (SpecialCase) {
                    if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
                    nRows = testRows;
                } else {
                    nRows = RUserCap;
                    if (nRows > testRows) {nRows = testRows;}
                }
            }
            
            XPtr<funcPtr> xpFun1 = putFunPtrInXPtr(mainFun1);
            funcPtr myFun1 = *xpFun1;
            
            if (SpecialCase) {
                NumericMatrix matRes;
                double testVal;
                bool test;
                std::vector<int> indexMatch;
                std::vector<double> matchVals;
                indexMatch.reserve(nRows);
                matchVals.reserve(nRows);
                
                if (IsComb) {
                    matRes = ComboGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes);
                } else {
                    matRes = PermuteGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes);
                }
                
                XPtr<compPtr> xpComp = putCompPtrInXPtr(compFun);
                compPtr myComp = *xpComp;
                
                for (i = 0; i < nRows; i++) {
                    for (j = 0; j < m; j++) {rowVec[j] = matRes(i, j);}
                    testVal = myFun1(rowVec, m);
                    test = myComp(testVal, myLim);
                    if (test) {
                        indexMatch.push_back(i);
                        matchVals.push_back(testVal);
                    }
                }
                
                int numCols = m;
                if (keepRes) {numCols++;}
                testRows = indexMatch.size();
                if (testRows > RUserCap) {nRows = RUserCap;} else {nRows = testRows;}
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
                                                               mainFun1, compFun, myLim, nRows, IsComb, true);
                int nRows = matRes.nrow();
                
                for (i = 0; i < nRows; i++) {
                    for (j = 0; j < m; j++) {vNum[j] = matRes(i, j);}
                    matRes(i, m) = myFun1(vNum, m);
                }
                
                return matRes;
            } else {
                if (IsInteger) {
                    return CombinatoricsConstraints<IntegerMatrix>(n, m, vNum, IsRepetition,
                                                                   mainFun1, compFun, myLim, nRows, IsComb, false);
                } else {
                    return CombinatoricsConstraints<NumericMatrix>(n, m, vNum, IsRepetition,
                                                                   mainFun1, compFun, myLim, nRows, IsComb, false);
                }
            }
        } else {
            if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
            nRows = testRows;
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
                    matRes = ComboGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, true);
                } else {
                    matRes = PermuteGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, true);
                }

                for (i = 0; i < nRows; i++) {
                    for (j = 0; j < m; j++) {rowVec[j] = matRes(i, j);}
                    matRes(i, m) = myFun2(rowVec, m);
                }
                
                return matRes;
            } else {
                if (IsComb) {
                    if (IsInteger) {
                        return ComboGeneral<IntegerMatrix>(n, m, vInt, IsRepetition, nRows, false);
                    } else if (IsFactor) {
                        return ComboFactor(n, m, vFactor, IsRepetition, nRows);
                    } else {
                        return ComboGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, false);
                    }
                } else {
                    if (IsInteger) {
                        return PermuteGeneral<IntegerMatrix>(n, m, vNum, IsRepetition, nRows, false);
                    } else if (IsFactor) {
                        return PermuteFactor(n, m, vFactor, IsRepetition, nRows);
                    } else {
                        return PermuteGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, false);
                    }
                }
            }
        }
    }
}