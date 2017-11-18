#include <Rcpp.h>
using namespace Rcpp;

std::vector<double> SectionLength(std::vector<int> v, std::vector<int> myRep) {
    unsigned long int n = v.size(), numRow = 1;
    unsigned long int i, j, k, mySum, n1 = n-1;
    double res;
    std::vector<int>::iterator it;
    std::vector<int> temp(n);
    for (it = v.begin(); it < v.end(); it++) {numRow *= (*it+1);}
    std::vector<double> myLengths(numRow);
    
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
    unsigned long int uN = n, count, strt, ind, divTest;
    unsigned long int i, j, k, m = 1, r, numCols;
    
    numCols = std::accumulate(Reps.begin(), Reps.end(), 0);
    unsigned long int lastCol = numCols - 1, lastSum = 0;
    
    int myInt;
    std::vector<int> repLen(n+1, 1);
    
    for (myInt = (n-1); myInt >= 0; myInt--) {
        repLen[myInt] = repLen[myInt+1]*(Reps[myInt]+1);
    }
    
    std::vector<double> groupLen = SectionLength(Reps, repLen);
    unsigned long int gLen = groupLen.size();
    
    if (groupLen[gLen - 1] > 2147483647) {
        stop("The number of rows cannot exceed 2^31 - 1.");
    }
    
    unsigned long int uRowN = groupLen[gLen - 1];
    TypeRcpp permuteMatrix(uRowN, numCols);
    
    std::vector<unsigned long int> vecLast(uRowN);
    
    std::vector<int> prevPerms(2*uRowN/5);
    std::vector<int> nextPerms(4*uRowN/5);
    typename std::vector<stdType>::iterator it, vBeg, vEnd;
    vBeg = v.begin(); vEnd = v.end();
    numCols--;
    nextPerms[0] = gLen - 1;
    
    for (i = 1; i < uN; i++) {
        for (j = 0; j < Reps[i]; j++) {
            lastSum += i;
        }
    }
    
    for (i = 0; i < uRowN; i++) {vecLast[i] = lastSum;}
    
    for (i = 0; i < (lastCol - 1); i++) {
        for (j = 0; j < m; j++) {prevPerms[j] = nextPerms[j];}
        strt = count = m = r = 0;
        while (count < uRowN) {
            ind = prevPerms[r];
            j = 0;
            for (it = vBeg; it < vEnd; it++) {
                divTest = (ind % repLen[j])/repLen[j+1];
                if (divTest > 0) {
                    ind -= repLen[j+1];
                    count += groupLen[ind];
                    for (k = strt; k < count; k++) {
                        permuteMatrix(k, i) = *it;
                        vecLast[k] -= j;
                    }
                    strt = count;
                    nextPerms[m] = ind;
                    ind += repLen[j+1];
                    m++;
                }
                j++;
            }
            r++;
        }
    }
    
    for (i = (lastCol-1); i < lastCol; i++) {
        strt = count = m = r = 0;
        while (count < uRowN) {
            ind = nextPerms[r];
            j = 0;
            for (it = vBeg; it < vEnd; it++) {
                divTest = (ind % repLen[j])/repLen[j+1];
                if (divTest > 0) {
                    ind -= repLen[j+1];
                    count += groupLen[ind];
                    for (k = strt; k < count; k++) {
                        permuteMatrix(k, i) = *it;
                        vecLast[k] -= j;
                    }
                    strt = count;
                    ind += repLen[j+1];
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
SEXP PermsRepLenRcpp(SEXP Rv, SEXP Rm, SEXP RIsFactor) {
    
    int i, j, m1, m2, lenR, lenV;
    std::vector<int> myReps;
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
            stop("each element in repLen must be positive");
        }
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