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

double NumPermsWithRep(std::vector<double> v) {
    List myRle = rleCpp(v);
    unsigned long int n = v.size(), myMax;
    std::vector<unsigned long int> myLens = myRle[0], myUnis = myRle[2];
    std::sort(myLens.begin(), myLens.end(),
              std::greater<unsigned long int>());
    
    myMax = myLens[0];
    unsigned long int i, j, numUni = myUnis[0];
    double result = 1;
    
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
    
    double testVal;
    int count = 0, vSize = v.size(), numCols;
    
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
        comparisonFunTwo = *xpCompOne;
    } else {
        std::sort(v.begin(), v.end());
        if (myComparison == "==") {
            XPtr<compPtr> xpCompThree = putCompPtrInXPtr("<=");
            comparisonFunTwo = *xpCompThree;
        } else {
            comparisonFunTwo = *xpCompOne;
        }
    }
    
    if (xtraCol) {numCols = r+1;} else {numCols = r;}
    TypeRcpp combinatoricsMatrix(rowNum, numCols);
    
    std::vector<double> z, testVec(r), zPerm(r);
    bool t_1, t_2, t = true, keepGoing = true;
    int r1 = r - 1, r2 = r - 2, maxZ = n - 1;
    int numPerms, k, i, j;
    
    if (repetition) {
        
        v.erase(std::unique(v.begin(), v.end()), v.end());
        vSize = v.size();
        z.assign(r, 0);
        maxZ = vSize - 1;
        
        while (keepGoing) {
            
            t_2 = true;
            for (i = 0; i < r; i++) {testVec[i] = v[z[i]];}
            testVal = constraintFun(testVec, r);
            t = comparisonFunTwo(testVal, lim);
            
            while (t && t_2 && keepGoing) {
                
                testVal = constraintFun(testVec, r);
                t_1 = comparisonFunOne(testVal, lim);
                
                if (t_1) {
                    if (isComb) {
                        for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = v[z[k]];}
                        count++;
                    } else {
                        zPerm = z;
                        numPerms = NumPermsWithRep(zPerm);
                        for (i=0; i < numPerms; i++) {
                            for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = v[zPerm[k]];}
                            std::next_permutation(zPerm.begin(), zPerm.end());
                            count++;
                            if (count >= rowNum) {break;}
                        }
                    }
                    keepGoing = (count < rowNum);
                }
                
                t_2 = (z[r1] != maxZ);
                
                if (t_2) {
                    z[r1]++;
                    testVec[r1] = v[z[r1]];
                    testVal = constraintFun(testVec, r);
                    t = comparisonFunTwo(testVal, lim);
                }
            }
            
            if (keepGoing) {
                for (i = r2; i >= 0; i--) {
                    if (z[i] != maxZ) {
                        z[i]++;
                        testVec[i] = v[z[i]];
                        for (k = (i+1); k < r; k++) {
                            z[k] = z[k-1];
                            testVec[k] = v[z[k]];
                        }
                        testVal = constraintFun(testVec, r);
                        t = comparisonFunTwo(testVal, lim);
                        if (t) {break;}
                    }
                }
                if (!t) {keepGoing = false;}
            }
        }
    } else {
        
        for (i = 0; i < r; i++) {z.push_back(i);}
        int indexRows = 0;
        IntegerMatrix indexMatrix;
        if (!isComb) {
            indexRows = (int)NumPermsNoRep(r, r-1);
            indexMatrix = MakeIndexHeaps(indexRows, r);
        }

        while (keepGoing) {
            
            t_2 = true;
            for (i = 0; i < r; i++) {testVec[i] = v[z[i]];}
            testVal = constraintFun(testVec, r);
            t = comparisonFunTwo(testVal, lim);
            
            while (t && t_2 && keepGoing) {
                
                testVal = constraintFun(testVec, r);
                t_1 = comparisonFunOne(testVal, lim);
                
                if (t_1) {
                    if (isComb) {
                        for (k=0; k < r; k++) {combinatoricsMatrix(count, k) = v[z[k]];}
                        count++;
                    } else {
                        for (j = 0; j < indexRows; j++) {
                            for (k = 0; k < r; k++) {
                                combinatoricsMatrix(count, k) = v[z[indexMatrix(j, k)]];
                            }
                            count++;
                            if (count >= rowNum) {break;}
                        }
                    }
                    keepGoing = (count < rowNum);
                }
                
                t_2 = (z[r1] != maxZ);
                
                if (t_2) {
                    z[r1]++;
                    testVec[r1] = v[z[r1]];
                    testVal = constraintFun(testVec, r);
                    t = comparisonFunTwo(testVal, lim);
                }
            }
            
            if (keepGoing) {
                for (i = r2; i >= 0; i--) {
                    if (z[i] != (n - r + i)) {
                        z[i]++;
                        testVec[i] = v[z[i]];
                        for (k = (i+1); k < r; k++) {
                            z[k] = z[k-1]+1;
                            testVec[k] = v[z[k]];
                        }
                        testVal = constraintFun(testVec, r);
                        t = comparisonFunTwo(testVal, lim);
                        if (t) {break;}
                    }
                }
                if (!t) {keepGoing = false;}
            }
        }
    }
    return(SubMat(combinatoricsMatrix, count));
}

template <typename TypeRcpp, typename stdType>
TypeRcpp ComboGeneral(int n, int r, std::vector<stdType> v,
                      bool repetition, int rowNum, bool xtraCol) {
    std::sort(v.begin(), v.end());
    std::vector<int> z;
    int r1 = r - 1, r2 = r - 2;
    int k, i, numIter, vSize;
    int numCols, maxZ, count = 0;
    if (xtraCol) {numCols  = r + 1;} else {numCols = r;}
    TypeRcpp combinationMatrix(rowNum, numCols);
    
    if (repetition) {
        v.erase(std::unique(v.begin(), v.end()), v.end());
        vSize = v.size();
        z.assign(r, 0);
        maxZ = vSize - 1;

        while (count < rowNum) {
            numIter = vSize - z[r1];
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
        while (count < rowNum) {
            numIter = n - z[r1];
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
        unsigned long int groupLen = 1, repLen = 1, myCol;
        typename std::vector<stdType>::iterator m, vBeg, vEnd;
        vBeg = v.begin(); vEnd = v.end();
        for (i = uR; i > 0; i--) {
            myCol = i-1;
            groupLen *= uN;
            chunk = 0;
            for (k = 0; k < uRowN; k += groupLen) {
                for (m = vBeg; m < vEnd; m++) {
                    for (j = 0; j < repLen; j++) {
                        permuteMatrix(chunk + j, myCol) = *m;
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

template <typename TypeRcpp, typename stdType>
TypeRcpp MultisetPermutation(int n, std::vector<stdType> v,
                             std::vector<int> Reps, int numRows) {
    
    int i, j = 0, numCols = 0;
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
    std::vector<stdType> permVec;
    typename std::vector<stdType>::iterator myIter, permEnd;
    TypeRcpp permuteMatrix(numRows, numCols);
    
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
    
    return permuteMatrix;
}

// [[Rcpp::export]]
SEXP CombinatoricsRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, 
                       SEXP f1, SEXP f2, SEXP lim, SEXP numRow,
                       SEXP RIsComb, SEXP RIsFactor,
                       SEXP RKeepRes, SEXP RFreqs) {
    
    int n, m = 0, i, j, m1, m2;
    int lenFreqs = 0, nRows = 0;
    double testRows, seqEnd;
    bool IsRepetition, IsInteger;
    bool keepRes, IsComb, IsFactor;
    bool SpecialCase, IsCharacter;
    bool IsConstrained, IsMultiset = false;
    
    std::vector<double> vNum, freqsExpanded;
    IntegerVector vFactor;
    std::vector<int> vInt, myReps;
    std::vector<std::string > vStr;
    
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
    
    IsComb = as<bool>(RIsComb);
    
    if (Rf_isNull(Rm)) {
        if (Rf_isNull(RFreqs)) {
            m = Rf_length(Rv);
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
    
    if (!IsMultiset) {if (m < 1) {stop("m must be positive");}}
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
        n = vNum.size();
    }
     
    if (IsInteger) {vInt.assign(vNum.begin(), vNum.end());}
    if (IsFactor) {IsCharacter = IsInteger = false;}

    if (IsMultiset) {
        if (n != lenFreqs) {stop("the length of freqs must equal the length of v");}
        testRows = NumPermsWithRep(freqsExpanded);
    } else {
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
    }

    if (IsCharacter) {
        if (testRows > 2147483647) {stop("The number of rows cannot exceed 2^31 - 1.");}
        nRows = testRows;
        if (IsMultiset) {
            return MultisetPermutation<CharacterMatrix>(n, vStr, myReps, nRows);
        } else {
            if (IsComb){
                return ComboGeneral<CharacterMatrix>(n, m, vStr, IsRepetition, nRows, false);
            } else {
                return PermuteGeneral<CharacterMatrix>(n, m, vStr, IsRepetition, nRows, false);
            }
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
                } else if (IsMultiset) {
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
                    myLim = as<double>(lim);
                    break;
                }
                default: {
                    stop("limitConstraints must be of type numeric or integer");
                }
            }

            if (NumericVector::is_na(myLim)) {stop("limitConstraints cannot be NA");}

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

            for (i = (n-1); i >= 0; i--) {
                if (NumericVector::is_na(vNum[i])) {
                    vNum.erase(vNum.begin() + i);
                }
            }

            j = vNum.size();

            if (j < n) {
                n = j;
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
                bool Success;
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
                    Success = myComp(testVal, myLim);
                    if (Success) {
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

            if (Rf_isNull(f1) || IsMultiset) {keepRes = false;}
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
                if (IsFactor) {
                    IntegerMatrix factorMat;
                    IntegerVector testFactor = as<IntegerVector>(Rv);
                    vFactor = seq(1, n);
                    vInt.assign(vFactor.begin(), vFactor.end());
                    CharacterVector myClass = testFactor.attr("class");
                    CharacterVector myLevels = testFactor.attr("levels");

                    if (IsMultiset) {
                        factorMat = MultisetPermutation<IntegerMatrix>(n, vInt, myReps, nRows);
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
                    if (IsMultiset) {
                        if (IsInteger) {
                            return MultisetPermutation<IntegerMatrix>(n, vInt, myReps, nRows);
                        } else {
                            return MultisetPermutation<NumericMatrix>(n, vNum, myReps, nRows);
                        }
                    } else {
                        if (IsComb) {
                            if (IsInteger) {
                                return ComboGeneral<IntegerMatrix>(n, m, vInt, IsRepetition, nRows, false);
                            } else {
                                return ComboGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, false);
                            }
                        } else {
                            if (IsInteger) {
                                return PermuteGeneral<IntegerMatrix>(n, m, vInt, IsRepetition, nRows, false);
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