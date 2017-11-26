#include <Rcpp.h>
using namespace Rcpp;

std::vector<unsigned long int> SectionLength(std::vector<int> v, std::vector<int> myRep) {
    unsigned long int n = v.size(), numRow = 1;
    unsigned long int i, j, k, mySum, n1 = n-1, res;
    std::vector<int>::iterator it;
    std::vector<int> temp(n);
    for (it = v.begin(); it < v.end(); it++) {numRow *= (*it+1);}
    std::vector<unsigned long int> myLengths(numRow);
    
    for (i = 0; i < numRow; i++) {
        temp[0] = i / myRep[1];
        mySum = temp[0];
        temp[n1] = (i % myRep[n1]);
        mySum += temp[n1];
        for (j = 1; j < n1; j++) {
            temp[j] = (i % myRep[j]) / myRep[j+1];
            mySum += temp[j];
        }
        std::sort(temp.begin(), temp.end(),
                  std::greater<unsigned long int>());
        res = 1;
        for (j = mySum; j > temp[0]; j--) {res *= j;}
        
        for (j = 1; j < n; j++) {
            if (temp[j] > 1) {
                for (k = 2; k <= temp[j]; k++) {
                    res /= k;
                }
            }
        }
        myLengths[i] = res;
    }
    
    return myLengths;
}

template <typename TypeRcpp, typename stdType>
TypeRcpp PermuteSpecificReps(int n, std::vector<stdType> v, std::vector<int> Reps) {
    unsigned long int count, indTemp, ind, divTest, uN = n;
    unsigned long int i, j, k, strt, numCols, m = 1;
    
    numCols = std::accumulate(Reps.begin(), Reps.end(), 0);
    unsigned long int lastSum = 0, pentUlt;
    
    int myInt;
    std::vector<int> repLen(n+1, 1);
    
    for (myInt = (n-1); myInt >= 0; myInt--) {
        repLen[myInt] = repLen[myInt+1]*(Reps[myInt]+1);
    }
    
    std::vector<unsigned long int> groupLen = SectionLength(Reps, repLen);
    unsigned long int gLen = groupLen.size();
    typename std::vector<stdType>::iterator it, vBeg, vEnd;
    vBeg = v.begin(); vEnd = v.end();
    
    unsigned long int uRowN = groupLen[gLen - 1];
    TypeRcpp permuteMatrix(uRowN, numCols);
    std::vector<unsigned short int> vecLast(uRowN), indexOne(uRowN), indexTwo(uRowN);
    std::vector<unsigned short int>::iterator uLit, uLEnd;
    numCols--;
    pentUlt = numCols - 1;
    
    for (i = 1; i < uN; i++) {
        for (j = 0; j < Reps[i]; j++) {
            lastSum += i;
        }
    }
    
    for (i = 0; i < uRowN; i++) {vecLast[i] = lastSum;}
    indexOne[0] = gLen - 1;
    
    for (i = 0; i < pentUlt; i++) {
        if (i % 2 == 0) {
            uLEnd = indexOne.begin() + m;
            strt = count = m = 0;
            for (uLit = indexOne.begin(); uLit < uLEnd; uLit++) {
                ind = *(uLit);
                j = 0;
                for (it = vBeg; it < vEnd; it++) {
                    divTest = (ind % repLen[j]);
                    if (divTest >= repLen[j+1]) {
                        indTemp = ind;
                        ind -= repLen[j+1];
                        count += groupLen[ind];
                        for (k = strt; k < count; k++) {
                            permuteMatrix(k, i) = *it;
                            vecLast[k] -= j;
                        }
                        strt = count;
                        indexTwo[m] = ind;
                        ind = indTemp;
                        m++;
                    }
                    j++;
                }
            }
        } else {
            uLEnd = indexTwo.begin() + m;
            strt = count = m = 0;
            for (uLit = indexTwo.begin(); uLit < uLEnd; uLit++) {
                ind = *(uLit);
                j = 0;
                for (it = vBeg; it < vEnd; it++) {
                    divTest = (ind % repLen[j]);
                    if (divTest >= repLen[j+1]) {
                        indTemp = ind;
                        ind -= repLen[j+1];
                        count += groupLen[ind];
                        for (k = strt; k < count; k++) {
                            permuteMatrix(k, i) = *it;
                            vecLast[k] -= j;
                        }
                        strt = count;
                        indexOne[m] = ind;
                        ind = indTemp;
                        m++;
                    }
                    j++;
                }
            }
        }
    }
    
    if (pentUlt % 2 == 0) {
        uLEnd = indexOne.begin() + m;
        strt = count = 0;
        for (uLit = indexOne.begin(); uLit < uLEnd; uLit++) {
            ind = *(uLit);
            j = 0;
            for (it = vBeg; it < vEnd; it++) {
                divTest = (ind % repLen[j]);
                if (divTest >= repLen[j+1]) {
                    indTemp = ind;
                    ind -= repLen[j+1];
                    count += groupLen[ind];
                    for (k = strt; k < count; k++) {
                        permuteMatrix(k, pentUlt) = *it;
                        vecLast[k] -= j;
                        permuteMatrix(k, numCols) = v[vecLast[k]];
                    }
                    strt = count;
                    ind = indTemp;
                }
                j++;
            }
        }
    } else {
        uLEnd = indexTwo.begin() + m;
        strt = count = m = 0;
        for (uLit = indexTwo.begin(); uLit < uLEnd; uLit++) {
            ind = *(uLit);
            j = 0;
            for (it = vBeg; it < vEnd; it++) {
                divTest = (ind % repLen[j]);
                if (divTest >= repLen[j+1]) {
                    indTemp = ind;
                    ind -= repLen[j+1];
                    count += groupLen[ind];
                    for (k = strt; k < count; k++) {
                        permuteMatrix(k, pentUlt) = *it;
                        vecLast[k] -= j;
                        permuteMatrix(k, numCols) = v[vecLast[k]];
                    }
                    strt = count;
                    ind = indTemp;
                }
                j++;
            }
        }
    }
    
    return(permuteMatrix);
}

double NumPermsWithRep(std::vector<int> v, double myMax) {
    std::sort(v.begin(), v.end(), std::greater<unsigned long int>());
    unsigned long int i, j, numUni = v.size();
    double result = 1;
    
    for (i = myMax; i > v[0]; i--) {result *= i;}
    
    if (numUni > 1) {
        for (i = 1; i < numUni; i++) {
            for (j = 2; j <= v[i]; j++) {
                result /= j;
            }
        }
    }
    
    return result;
}

// [[Rcpp::export]]
SEXP PermsRepLenRcpp(SEXP Rv, SEXP Rm, SEXP RIsFactor) {
    
    int i, m1, m2, lenR, lenV;
    std::vector<int> myReps;
    double seqEnd, rowTest = 0;
    bool IsCharacter, IsInteger, IsFactor;
    
    switch(TYPEOF(Rm)) {
        case REALSXP: {
            myReps = as<std::vector<int> >(Rm);
            break;
        }
        case INTSXP: {
            myReps = as<std::vector<int> >(Rm);
            break;
        }
        default: {
            stop("repLen must be of type numeric or integer");
        }
    }
    
    lenR = myReps.size();
    for (i = 0; i < lenR; i++) {
        if (myReps[i] < 1) {
            stop("each element in repLen must be a positive number");
        }
        rowTest += myReps[i];
    }
    
    rowTest = NumPermsWithRep(myReps, rowTest);
    if (rowTest > 2147483647) {
        stop("The number of rows cannot exceed 2^31 - 1.");
    }
    
    std::vector<double> vNum;
    IntegerVector vFactor;
    std::vector<int> vInt;
    std::vector<std::string > vStr;
    IsFactor = as<bool>(RIsFactor);

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
    
    if (IsCharacter) {
        vStr = as<std::vector<std::string > >(Rv);
        lenV = vStr.size();
        if (lenV == 1 && lenR == 1) {
            CharacterVector strVec(myReps[0], vStr[0]);
            return strVec;
        }
    } else {
        if (Rf_length(Rv) == 1) {
            seqEnd = as<double>(Rv);
            if (NumericVector::is_na(seqEnd)) {seqEnd = 1;}
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            IntegerVector vTemp = seq(m1, m2);
            IsInteger = true;
            vNum = as<std::vector<double> >(vTemp);
            
            if (lenR == 1) {
                if (IsFactor) {
                    IntegerVector facVec = as<IntegerVector>(Rv);
                    CharacterVector classVec = facVec.attr("class");
                    CharacterVector levelVec = facVec.attr("levels");
                    
                    IntegerVector facRet(myReps[0], 1);
                    facRet.attr("class") = classVec;
                    facRet.attr("levels") = levelVec;
                    
                    return facRet;
                } else {
                    if (IsInteger) {
                        IntegerVector intVec(myReps[0], seqEnd);
                        return intVec;
                    } else {
                        NumericVector numVec(myReps[0], seqEnd);
                        return numVec;
                    }
                }
            }
        } else {
            vNum = as<std::vector<double> >(Rv);
        }
        lenV = vNum.size();
    }

    if (lenR != lenV) {
        stop("The length of the source vector, v, "
                  "must equal the length of repLens.");
    }

    if (IsInteger) {vInt.assign(vNum.begin(), vNum.end());}
    if (IsFactor) {IsCharacter = IsInteger = false;}

    if (IsCharacter) {
        return PermuteSpecificReps<CharacterMatrix>(lenV, vStr, myReps);
    } else {
        if (IsInteger) {
            return PermuteSpecificReps<IntegerMatrix>(lenV, vInt, myReps);
        } else if (IsFactor) {
            IntegerVector testFactor = as<IntegerVector>(Rv);
            vFactor = seq(1, lenV);
            
            CharacterVector myClass = testFactor.attr("class");
            CharacterVector myLevels = testFactor.attr("levels");
            vInt.assign(vFactor.begin(), vFactor.end());
            
            IntegerMatrix factorMat;
            factorMat = PermuteSpecificReps<IntegerMatrix>(lenV, vInt, myReps);
            factorMat.attr("class") = myClass;
            factorMat.attr("levels") = myLevels;

            return factorMat;
        } else {
            return PermuteSpecificReps<NumericMatrix>(lenV, vNum, myReps);
        }
    }
}