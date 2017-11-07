#include <Rcpp.h>
#include <algorithm>
#include <string>
using namespace Rcpp;

List rleCpp(std::vector<double> x) {
    std::vector<unsigned long int> lengths;
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
    return List::create(
        _["lengths"] = lengths,
        _["values"] = values,
        _["uniques"] = i
    );
}

int GetNumPerms(std::vector<double> v) {
    List myRle = rleCpp(v);
    unsigned long int n = v.size(), myMax;
    std::vector<unsigned long int> myLens = myRle[0];
    std::sort(myLens.begin(), myLens.end(),
              std::greater<unsigned long int>());
    
    myMax = myLens[0];
    unsigned long int i, j, numUni = myRle[2];
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

double numPermutations(int n, int k) {
    double dblN = (double)n, result = 1;
    double i, m = dblN - (double)k;
    for (i = n; i > m; i--) {result *= i;}
    return result;
}

NumericMatrix SubMat(NumericMatrix m, int n) {
    int k = m.ncol();
    NumericMatrix subMatrix(n,k);
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

NumericMatrix CombinatoricsConstraints(int n, int r, std::vector<double> v,
                                       bool repetition, std::string myFun, 
                                       std::string myComparison, double lim,
                                       int rowNum, bool isComb) {
    // where myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max"
    // and myComparison is a comparison operator: "<", "<=", ">", or ">="
    
    double extremeV, testVal;
    int count = 0, s = v.size();
    unsigned int ur = r;
    
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
    
    NumericMatrix combinatoricsMatrix(rowNum, r);
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
                        numPerms = GetNumPerms(zPart);
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
        int indexRows = (int)numPermutations(r, r-1);
        IntegerMatrix indexMatrix = MakeIndexHeaps(indexRows, r);
        
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

NumericMatrix ComboNumeric(int n, int r, std::vector<double> v,
                                   bool repetition, int rowNum) {
    std::sort(v.begin(), v.end());
    std::vector<double> z;
    int count = 0, tSize = r - 1, s = v.size();
    int k, i, numIter, posMaxZ, newSize;
    double maxV = v[s-1];
    std::vector<double>::iterator it;
    std::vector<double> temp;
    bool keepGoing = true;
    NumericMatrix combinationMatrix(rowNum, r);
    
    if (repetition) {
        z.assign(r+n-1, v[0]);
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

CharacterMatrix ComboCharacter(int n, int r, std::vector<std::string > v,
                               bool repetition, int rowNum) {
    std::sort(v.begin(), v.end());
    std::vector<std::string > z;
    int count = 0, tSize = r - 1, s = v.size();
    int k, i, numIter, posMaxZ, newSize;
    std::string maxV = v[s-1];
    std::vector<std::string >::iterator it;
    std::vector<std::string > temp;
    bool keepGoing = true;
    CharacterMatrix combinationMatrix(rowNum, r);
    
    if (repetition) {
        z.assign(r+n-1, v[0]);
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

NumericMatrix PermuteNumeric(int n, int r, std::vector<double> v,
                             bool repetition, int rowNum) {
    unsigned long int uN = n, uR = r, uRowN = rowNum;
    unsigned long int i = 0, j, k, chunk;
    NumericMatrix permuteMatrix(uRowN, uR);
    
    if (repetition) {
        unsigned long int groupLen = 1, repLen = 1;
        std::vector<double>::iterator m, vBeg, vEnd;
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
        unsigned long int combRows = (int)nChooseK(n, r);
        NumericMatrix myCombs = ComboNumeric(n,r,v,false,combRows);
        int indexRows = (int)numPermutations(r, r-1);
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

CharacterMatrix PermuteCharacter(int n, int r, std::vector<std::string > v,
                                 bool repetition, int rowNum) {
    unsigned long int i = 0, j, k, chunk;
    unsigned long int uN = n, uR = r, uRowN = rowNum;
    CharacterMatrix permuteMatrix(uRowN, uR);
    
    if (repetition) {
        unsigned long int groupLen = 1, repLen = 1;
        std::vector<std::string >::iterator m, vBeg, vEnd;
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
        CharacterMatrix myCombs = ComboCharacter(n,r,v,false,combRows);
        int indexRows = (int)numPermutations(r, r-1);
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

// [[Rcpp::export]]
SEXP CombinatoricsRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, SEXP fun1,
                          SEXP fun2, SEXP lim, SEXP numRow, SEXP RIsComb) {
    int n, m, j, nRows;
    double testRows;
    bool IsRepetition, IsCharacter, IsConstrained, IsComb;
    
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
    std::vector<std::string > vStr;
    
    switch(TYPEOF(Rv)) {
        case REALSXP: {
            IsCharacter = false;
            break;
        }
        case INTSXP: {
            IsCharacter = false;
            break;
        }
        case STRSXP: {
            IsCharacter = true;
            break;
        }
        default: {
            stop("Only integers, numerical, and character classes are supported for v");   
        }
    }
    
    if (!Rf_isLogical(Rrepetition)) {stop("repetitions must be a logical value");}
    IsRepetition = as<bool >(Rrepetition);
    IsComb = as<bool>(RIsComb);
    
    if (IsCharacter) {
        vStr = as<std::vector<std::string > >(Rv);
        n = vStr.size();
    } else {
        if (Rf_length(Rv) == 1) {
            j = as<int>(Rv);
            if (j < m) {stop("v cannot be less than m");}
            IntegerVector vTemp = seq(1, j);
            vNum = as<std::vector<double> >(vTemp);
        } else {
            vNum = as<std::vector<double> >(Rv);    
        }
        n = vNum.size();
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
            testRows = numPermutations(n, m);
        }
    }
    
    if (IsCharacter) {
        if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
        nRows = testRows;
        if (IsComb){
            return ComboCharacter(n, m, vStr, IsRepetition, nRows);
        } else {
            return PermuteCharacter(n, m, vStr, IsRepetition, nRows);
        }
    } else {
        if (Rf_isNull(lim)) {
            IsConstrained = false;
        } else {
            if (Rf_isNull(fun1)) {
                IsConstrained = false;
            } else {
                if (!Rf_isString(fun1)) {stop("constraintFun must be passed as a character");}
                if (Rf_isNull(fun2)) {
                    IsConstrained = false;
                } else {
                    IsConstrained = true;
                    if (!Rf_isString(fun2)) {stop("comparisonFun must be passed as a character");}
                }
            }
        }
        
        if (IsConstrained) {
            double myLim, testRows2;
            int myRows;

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
                testRows2 = 0;
            } else {
                switch(TYPEOF(numRow)) {
                    case REALSXP: {
                        testRows2 = as<double>(numRow);
                        break;
                    }
                    case INTSXP: {
                        testRows2 = as<double>(numRow);
                        break;
                    }
                    default: {
                        stop("rowCap must be of type numeric or integer");
                    }
                }
            }

            if (testRows2 == 0) {
                if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
                myRows = testRows;
            } else if (testRows2 < 0) {
                stop("The number of rows must be positive");
            } else if (testRows2 > 2147483647) {
                stop("The number of rows cannot exceed 2^31 - 1.");
            } else {
                myRows = testRows2;
                if (myRows > testRows) {myRows = testRows;}
            }
            
            std::string mainFun = as<std::string >(fun1);
            if (mainFun != "prod" && mainFun != "sum" && mainFun != "mean"
                    && mainFun != "max" && mainFun != "min") {
                stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
            }

            std::string compFun = as<std::string >(fun2);
            if (compFun != "<" && compFun != "<=" && compFun != ">"
                    && compFun != ">=" && compFun != "==") {
                stop("comparisonFun must be one of the following: >, >=, <, <=, or ==");
            }

            return CombinatoricsConstraints(n, m, vNum, IsRepetition,
                                            mainFun, compFun, myLim, myRows, IsComb);
        } else {
            if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
            nRows = testRows;
            if (IsComb) {
                return ComboNumeric(n, m, vNum, IsRepetition, nRows);
            } else {
                return PermuteNumeric(n, m, vNum, IsRepetition, nRows);
            }
        }
    }
}