#include <Rcpp.h>
#include <algorithm>
#include <string>
using namespace Rcpp;

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

// This function applys a constraint function to a vector v with respect
// to a constraint value "lim". The main idea is that combinations are added
// successively, until a particular combination exceeds the given constraint
// value for a given constraint function. After this point, we can safely skip
// several combinations knowing that they will exceed the given constraint value.

NumericMatrix ComboConstraints(int n, int r, std::vector<double> v, bool repetition, 
                               std::string myFun, std::string myComparison, double lim, int rowNum) {
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
    
    NumericMatrix combinationMatrix(rowNum, r);
    std::vector<double> z;
    bool t_1, t_2, t = true, keepGoing = true;
    std::vector<double>::iterator j;
    std::vector<double> temp;
    int posMaxZ, tSize = r - 1, newSize, rMinus, k, zSize;
    
    if (repetition) {
        z.assign(r+n-1, v[0]);
        std::copy(v.begin()+1, v.end(), z.begin()+r);
        
        while (keepGoing) {
            t_2 = z.size() >= ur;
            testVal = constraintFun(z, r);
            t = comparisonFunTwo(testVal, lim);
            while (t && t_2 && keepGoing) {
                testVal = constraintFun(z, r);
                t_1 = comparisonFunOne(testVal, lim);
                if (t_1) {
                    for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
                    count++;
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
                        for (j = z.begin(); j <= z.begin() + posMaxZ; j++) {temp.push_back(*j);}
        
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
                            for (j = z.begin(); j < z.begin() + rMinus; j++) {temp.push_back(*j);}
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
        while (keepGoing) {
            t_2 = z.size() >= ur;
            testVal = constraintFun(z, r);
            t = comparisonFunTwo(testVal, lim);
            while (t && t_2 && keepGoing) {
                testVal = constraintFun(z, r);
                t_1 = comparisonFunOne(testVal, lim);
                if (t_1) {
                    for (k=0; k < r; k++) {combinationMatrix(count, k) = z[k];}
                    count++;
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
                        for (j = z.begin(); j <= z.begin() + posMaxZ; j++) {temp.push_back(*j);}
                        for (j = v.begin() + (k+2); j < v.end(); j++) {temp.push_back(*j);}
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
                            for (j = z.begin(); j <= z.begin() + rMinus; j++) {temp.push_back(*j);}
                            for (j = v.begin() + (k+2); j < v.end(); j++) {temp.push_back(*j);}
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
    return(SubMat(combinationMatrix, count));
}

NumericMatrix ComboStandardNumeric(int n, int r, std::vector<double> v,
                                   bool repetition, int rowNum) {
    std::sort(v.begin(), v.end());
    std::vector<double> z;
    int count = 0, tSize = r - 1, s = v.size();
    int k, i, numIter, posMaxZ, newSize;
    double maxV = v[s-1];
    std::vector<double>::iterator j;
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
                for (j = z.begin(); j <= z.begin() + posMaxZ; j++) {temp.push_back(*j);}
    
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
                for (j = z.begin(); j <= z.begin() + posMaxZ; j++) {temp.push_back(*j);}
                for (j = v.begin() + (k+2); j < v.end(); j++) {temp.push_back(*j);}
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
    std::vector<std::string >::iterator j;
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
                for (j = z.begin(); j <= z.begin() + posMaxZ; j++) {temp.push_back(*j);}
    
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
                for (j = z.begin(); j <= z.begin() + posMaxZ; j++) {temp.push_back(*j);}
                for (j = v.begin() + (k+2); j < v.end(); j++) {temp.push_back(*j);}
                z = temp;
            } else {
                keepGoing = false;
            }
        }
    }
    return(combinationMatrix);
}

// [[Rcpp::export]]
SEXP ComboRcpp(SEXP Rv, SEXP Rm, SEXP Rrep, SEXP fun1,
                          SEXP fun2, SEXP lim, SEXP numRow) {
    int n, m, j, nRows;
    double testRows;
    bool rep, IsCharacter, IsConstrained;
    
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
    
    if (!Rf_isLogical(Rrep)) {stop("repetitions must be a logical value");}
    rep = as<bool >(Rrep);
    
    if (IsCharacter) {
        std::vector<std::string > vStr = as<std::vector<std::string > >(Rv);
        n = vStr.size();
        if (m > n) {stop("m must be less than or equal to the length of v");}
        if (rep) {
            testRows = GetRowNum(n, m);    
        } else {
            testRows = nChooseK(n, m);
        }
        if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
        nRows = testRows;
        return ComboCharacter(n, m, vStr, rep, nRows);
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
        if (m > n) {stop("m must be less than or equal to the length of v");}
        
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
        
        if (rep) {
            testRows = GetRowNum(n, m);    
        } else {
            testRows = nChooseK(n, m);
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
            
            return ComboConstraints(n, m, vNum, rep, mainFun, compFun, myLim, myRows);
        } else {
            if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
            nRows = testRows;
            return ComboStandardNumeric(n, m, vNum, rep, nRows);
        }
    }
}