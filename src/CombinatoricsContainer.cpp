#include <Combinations.h>
#include <Permutations.h>
#include <Combinatorics.h>
#include <NthResult.h>
using namespace Rcpp;

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
    
    // myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max"
    // and myComparison is a comparison operator: "<", "<=", ">", or ">="
    
    double testVal;
    int count = 0, vSize = v.size(), numCols;
    int testRows, i, j, numRows2 = 0;
    
    XPtr<funcPtr> xpFun = putFunPtrInXPtr(myFun);
    funcPtr constraintFun = *xpFun;
    
    XPtr<compPtr> xpCompOne = putCompPtrInXPtr(myComparison);
    compPtr comparisonFunOne = *xpCompOne;
    
    XPtr<compPtr> xpCompTwo = xpCompOne;
    compPtr comparisonFunTwo;

    if (myComparison == ">" || myComparison == ">=") {
        if (isMult) {
            for (i = 0; i < (n - 1); i++) {
                for (j = (i + 1); j < n; j++) {
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
        if (testRows == numRows)
            executeCap = false;
    }
    
    numCols = xtraCol ? (r + 1) : r;
    typeRcpp combinatoricsMatrix(numRows, numCols);
    
    std::vector<int> z, zCheck, zPerm(r);
    std::vector<double> testVec(r);
    bool t_1, t_2, t = true, keepGoing = true;
    int r1 = r - 1, r2 = r - 2, k = 0; 
    int numIter, maxZ = n - 1;
    
    if (isMult) {
        int zExpSize = 0;
        std::vector<int> eachRowCount;
        if (!isComb)
            eachRowCount.reserve(numRows);
        
        std::vector<int> zExpand, zIndex, zGroup(r);
        
        for (i = 0; i < n; i++)
            zExpSize += Reps[i];
        
        zIndex.reserve(n);
        zExpand.reserve(zExpSize);
        for (i = 0; i < n; i++) {
            zIndex.push_back(k);
            for (j = 0; j < Reps[i]; j++, k++)
                zExpand.push_back(i);
        }
        
        z.reserve(r);
        for (i = 0; i < r; i++)
            z.push_back(zExpand[i]);
        
        while (keepGoing) {
            
            t_2 = true;
            for (i = 0; i < r; i++)
                testVec[i] = v[zExpand[zIndex[z[i]]]];
            
            testVal = constraintFun(testVec);
            t = comparisonFunTwo(testVal, lim);
            
            while (t && t_2 && keepGoing) {
                
                testVal = constraintFun(testVec);
                t_1 = comparisonFunOne(testVal, lim);
                
                if (t_1) {
                    for (k = 0; k < r; k++)
                        combinatoricsMatrix(count, k) = v[zExpand[zIndex[z[k]]]];
                    
                    count++;
                    
                    if (!isComb) {
                        for (k = 0; k < r; k++)
                            zPerm[k] = zExpand[zIndex[z[k]]];
                        
                        numIter = (int) NumPermsWithRep(zPerm);
                        numRows2 += numIter;
                        eachRowCount.push_back(numIter);
                    }
                    
                    if (executeCap)
                        keepGoing = (count < testRows);    
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
            int count2 = 0, segment = 0;
            bool bTooMany = false;
            
            for (i = 0; i < count; i++) {
                // populate permVec with next combination
                for (j = 0; j < r; j++)
                    testVec[j] = combinatoricsMatrix(i, j);
                
                segment = eachRowCount[i];
                if (executeCap) {
                    if ((eachRowCount[i] + count2) > testRows) {
                        segment = testRows - count2;
                        bTooMany = true;
                    }
                }
                
                // populate permuteMatrix with all permutations
                // of a particular combination
                for (j = 0; j < segment; j++, count2++) {
                    for (k = 0; k < r; k++)
                        permuteMatrix(count2, k) = testVec[k];

                    std::next_permutation(testVec.begin(), testVec.end());
                }
                if (bTooMany) {break;}
            }
            
            return Combinations::SubMat(permuteMatrix, count2);
        }
        
    } else if (repetition) {
        
        v.erase(std::unique(v.begin(), v.end()), v.end());
        vSize = v.size();
        z.assign(r, 0);
        maxZ = vSize - 1;
        
        while (keepGoing) {
            
            t_2 = true;
            for (i = 0; i < r; i++)
                testVec[i] = v[z[i]];
            
            testVal = constraintFun(testVec);
            t = comparisonFunTwo(testVal, lim);
            
            while (t && t_2 && keepGoing) {
                
                testVal = constraintFun(testVec);
                t_1 = comparisonFunOne(testVal, lim);
                
                if (t_1) {
                    if (isComb) {
                        for (k = 0; k < r; k++)
                            combinatoricsMatrix(count, k) = v[z[k]];
                        
                        count++;
                    } else {
                        zPerm = z;
                        numIter = (int) NumPermsWithRep(zPerm);
                        if ((numIter + count) > numRows)
                            numIter = numRows - count;
                        
                        for (i = 0; i < numIter; i++, count++) {
                            for (k = 0; k < r; k++)
                                combinatoricsMatrix(count, k) = v[zPerm[k]];
                            
                            std::next_permutation(zPerm.begin(), zPerm.end());
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
        
        for (i = 0; i < r; i++)
            z.push_back(i);
        
        int indexRows, nMinusR = (n - r), myRow;
        indexRows = isComb ? 0 : (int) NumPermsNoRep(r, r1);
        uint8_t *indexMatrix = new uint8_t[indexRows * r];
        
        if (!isComb) {
            indexRows = (int) NumPermsNoRep(r, r1);
            std::vector<uint8_t> indexVec(r);
            std::iota(indexVec.begin(), indexVec.end(), 0);
            
            for (i = 0, myRow = 0; i < indexRows; i++, myRow += r) {
                for (j = 0; j < r; j++)
                    indexMatrix[myRow + j] = indexVec[j];
                
                std::next_permutation(indexVec.begin(), indexVec.end());
            }
        }

        while (keepGoing) {

            t_2 = true;
            for (i = 0; i < r; i++)
                testVec[i] = v[z[i]];

            testVal = constraintFun(testVec);
            t = comparisonFunTwo(testVal, lim);

            while (t && t_2 && keepGoing) {

                testVal = constraintFun(testVec);
                t_1 = comparisonFunOne(testVal, lim);

                if (t_1) {
                    if (isComb) {
                        for (k=0; k < r; k++)
                            combinatoricsMatrix(count, k) = v[z[k]];

                        count++;
                    } else {
                        if (indexRows + count > numRows)
                            indexRows = numRows - count;

                        for (j = 0, myRow = 0; j < indexRows; j++, count++, myRow += r)
                            for (k = 0; k < r; k++)
                                combinatoricsMatrix(count, k) = v[z[indexMatrix[myRow + k]]];
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
        
        delete[] indexMatrix;
    }
       
    return(Combinations::SubMat(combinatoricsMatrix, count));
}

// [[Rcpp::export]]
SEXP CombinatoricsRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, SEXP RFreqs,
                       SEXP Rlow, SEXP Rhigh, SEXP f1, SEXP f2,
                       SEXP lim, bool IsComb, SEXP RKeepRes, 
                       bool IsFactor, bool IsCount) {
    
    int n, m = 0, m1, m2;
    int lenFreqs = 0, nRows = 0;
    double lower = 0, upper = 0, computedRows, seqEnd;
    bool IsRepetition, SpecialCase, IsConstrained;
    bool IsMultiset, keepRes, IsInteger, IsCharacter;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqsExpanded;
    CharacterVector vStr;
    
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
        
        lenFreqs = (int) myReps.size();
        for (int i = 0; i < lenFreqs; i++) {
            if (myReps[i] < 1) 
                stop("each element in freqs must be a positive number");
            
            for (int j = 0; j < myReps[i]; j++)
                freqsExpanded.push_back(i);
        }
    }
    
    if (Rf_isNull(Rm)) {
        if (IsMultiset) {
            m = freqsExpanded.size();
        } else {
            stop("m and freqs cannot both be NULL");
        }
    } else {
        if (Rf_length(Rm) > 1)
            stop("length of m must be 1");
        
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
    
    if (m < 1)
        stop("m must be positive");
    
    std::vector<double> rowVec(m);
    
    if (!Rf_isLogical(Rrepetition))
        stop("repetitions must be a logical value");
    
    IsRepetition = as<bool>(Rrepetition);
    
    if (!Rf_isLogical(RKeepRes))
        stop("keepResults must be a logical value");
    
    keepRes = as<bool>(RKeepRes);
    
    if (IsCharacter) {
        vStr = as<CharacterVector>(Rv);
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
        
        for (int i = (vNum.size() - 1); i >= 0; i--)
            if (NumericVector::is_na(vNum[i]))
                vNum.erase(vNum.begin() + i);
            
        n = vNum.size();
    }
    
    if (IsInteger)
        vInt.assign(vNum.begin(), vNum.end());
    
    if (IsFactor)
        IsCharacter = IsInteger = false;
    
    if (Rf_isNull(lim)) {
        IsConstrained = false;
    } else {
        if (Rf_isNull(f1)) {
            IsConstrained = false;
        } else {
            if (!Rf_isString(f1))
                stop("constraintFun must be passed as a character");
            
            if (Rf_isNull(f2)) {
                IsConstrained = false;
            } else if (IsFactor) {
                IsConstrained = false;
            } else {
                IsConstrained = true;
                if (!Rf_isString(f2))
                    stop("comparisonFun must be passed as a character");
            }
        }
    }
    
    bool bLower = false;
    bool bUpper = false;
    nRows = 0;
    
    if (!IsCount) {
        if (Rf_isNull(Rlow)) {
            lower = 0;
        } else { 
            bLower = true;
            switch(TYPEOF(Rlow)) {
                case REALSXP: {
                    lower = as<double>(Rlow);
                    break;
                }
                case INTSXP: {
                    lower = as<double>(Rlow);
                    break;
                }
                default: {
                    stop("bounds must be of type numeric or integer");
                }
            }
            
            lower--;
            if (lower < 0)
                stop("bounds must be positive");
        }
        
        if (Rf_isNull(Rhigh)) {
            upper = 0;
        } else {
            bUpper = true;
            switch(TYPEOF(Rhigh)) {
                case REALSXP: {
                    upper = as<double>(Rhigh);
                    break;
                }
                case INTSXP: {
                    upper = as<double>(Rhigh);
                    break;
                }
                default: {
                    stop("bounds must be of type numeric or integer");
                }
            }
        
            if (upper < 0)
                stop("bounds must be positive");
        }
    }
    
    if (IsMultiset) {
        if (n != lenFreqs)
            stop("the length of freqs must equal the length of v");
        
        if (IsComb || IsConstrained) {
            if (m > (int) freqsExpanded.size())
                m = freqsExpanded.size();
            
            computedRows = MultisetCombRowNum(n, m, myReps);
        } else {
            if (Rf_isNull(Rm)) {
                computedRows = NumPermsWithRep(freqsExpanded);
            } else if (m == (int) freqsExpanded.size()) {
                computedRows = NumPermsWithRep(freqsExpanded);
            } else {
                if (m > (int) freqsExpanded.size())
                    m = freqsExpanded.size();
                
                computedRows = MultisetPermRowNum(n, m, myReps);
            }
        }
    } else {
        if (IsRepetition) {
            if (IsComb)
                computedRows = NumCombsWithRep(n, m);
            else
                computedRows = std::pow((double) n, (double) m);
        } else {
            if (m > n)
                stop("m must be less than or equal to the length of v");
            
            if (IsComb)
                computedRows = nChooseK(n, m);
            else
                computedRows = NumPermsNoRep(n, m);
        }
    }

    if (lower >= computedRows)
        stop("lower bound exceeds the maximum number of possible results");
    
    if (upper > computedRows)
        stop("upper bound exceeds the maximum number of possible results");
    
    if (IsCount)
        return wrap(computedRows);
    
    std::vector<int> startZ(m);
    bool permNonTrivial = false;
    
    if (bLower && lower > 0) {
        if (IsComb) {
            startZ = nthCombination(n, m, lower, IsRepetition, IsMultiset, myReps);
        } else {
            permNonTrivial = true;
            startZ = nthPermutation(n, m, lower, IsRepetition, IsMultiset, myReps);
            
            if (IsMultiset) {
                
                for (std::size_t j = 0; j < startZ.size(); j++) {
                    for (std::size_t i = 0; i < freqsExpanded.size(); i++) {
                        if (freqsExpanded[i] == startZ[j]) {
                            freqsExpanded.erase(freqsExpanded.begin() + i);
                            break;
                        }
                    }
                }
                
                for (std::size_t i = 0; i < freqsExpanded.size(); i++)
                    startZ.push_back(freqsExpanded[i]);
                
            } else if (!IsRepetition) {
                
                if (m < n) {
                    for (int i = 0; i < n; i++) {
                        bool bExist = false;
                        for (std::size_t j = 0; j < startZ.size(); j++) {
                            if (startZ[j] == i) {
                                bExist = true;
                                break;
                            }
                        }
                        if (!bExist)
                            startZ.push_back(i);
                    }
                }
            }
        }
    } else {
        if (IsComb) {
            if (IsMultiset)
                startZ.assign(freqsExpanded.begin(), freqsExpanded.begin() + m);
            else if (IsRepetition)
                std::fill(startZ.begin(), startZ.end(), 0);
            else
                std::iota(startZ.begin(), startZ.end(), 0);
        } else {
            // Both repetition and non-repetition will
            // have an all zero vec, so there is no
            // need to have a condition for repetition
            if (IsMultiset) {
                startZ = freqsExpanded;
            } else {
                std::fill(startZ.begin(), startZ.end(), 0);
            }
        }
    }
    
    double userNumRows = 0;
    
    if (bLower && bUpper) {
        userNumRows = upper - lower;
    } else if (bUpper) {
        userNumRows = upper;
    } else if (bLower) {
        userNumRows = computedRows - lower;
    }
    
    if (userNumRows == 0) {
        if (computedRows > INT_MAX)
            stop("The number of rows cannot exceed 2^31 - 1.");
        
        userNumRows = nRows = computedRows;
    } else if (userNumRows < 0) {
        stop("The number of rows must be positive. Either the lowerBound "
              "exceeds the maximum number of possible results or the "
              "lowerBound is greater than the upperBound.");
    } else if (userNumRows > INT_MAX) {
        stop("The number of rows cannot exceed 2^31 - 1.");
    } else {
        nRows = userNumRows;
        if (nRows > computedRows)
            nRows = computedRows;
    }
    
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
        
        if (NumericVector::is_na(myLim))
            stop("limitConstraints cannot be NA");
        
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
            for (int i = 0; i < n; i++) {
                if (vNum[i] < 0) {
                    SpecialCase = true;
                    break;
                }
            }
        }
        
        XPtr<funcPtr> xpFun1 = putFunPtrInXPtr(mainFun1);
        funcPtr myFun1 = *xpFun1;
        
        if (SpecialCase) {
            if (computedRows > INT_MAX)
                stop("The number of rows cannot exceed 2^31 - 1.");
            
            nRows = (int) computedRows;
            NumericMatrix matRes;
            double testVal;
            bool Success;
            std::vector<int> indexMatch;
            std::vector<double> matchVals;
            indexMatch.reserve(nRows);
            matchVals.reserve(nRows);
            
            if (IsComb) {
                if (IsMultiset) {
                    matRes = Combinations::MultisetCombination<NumericMatrix>(n, m, vNum, myReps, nRows, keepRes, startZ);
                } else {
                    matRes = Combinations::ComboGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes, startZ);
                }
            } else {
                if (IsMultiset) {
                    matRes = Permutations::MultisetPermutation<NumericMatrix>(n, m, vNum, myReps, nRows, keepRes, startZ);
                    nRows = matRes.nrow();
                } else {
                    matRes = Permutations::PermuteGeneral<NumericMatrix>(n, m, vNum, IsRepetition, nRows, 
                                                                         keepRes, startZ, permNonTrivial);
                }
            }
            
            XPtr<compPtr> xpComp = putCompPtrInXPtr(compFun);
            compPtr myComp = *xpComp;
            
            for (int i = 0; i < nRows; i++) {
                for (int j = 0; j < m; j++)
                    rowVec[j] = matRes(i, j);
                
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
            nRows  = (computedRows > userNumRows) ? userNumRows : computedRows;
            NumericMatrix returnMatrix(nRows, numCols);
            
            for (int i = 0; i < nRows; i++) {
                for (int j = 0; j < m; j++)
                    returnMatrix(i,j) = matRes(indexMatch[i],j);

                if (keepRes)
                    returnMatrix(i,m) = matchVals[i];
            }
            
            return returnMatrix;
        }
        
        if (keepRes) {
            NumericMatrix matRes = CombinatoricsConstraints<NumericMatrix>(n, m, vNum, IsRepetition,
                                                                           mainFun1, compFun, myLim, nRows,
                                                                           IsComb, true, myReps, IsMultiset);
            nRows = matRes.nrow();
            
            for (int i = 0; i < nRows; i++) {
                for (int j = 0; j < m; j++)
                    rowVec[j] = matRes(i, j);
                
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
        if (IsComb) {
            return CombinationsRcpp(n, m, IsRepetition, vStr, nRows, vInt,
                                    vNum, IsMultiset, IsFactor, keepRes, startZ,
                                    IsCharacter, Rv, IsInteger, myReps, f1, f2);
        } else {
            return PermutationsRcpp(n, m, IsRepetition, vStr, nRows, vInt,
                                    vNum, IsMultiset, IsFactor, keepRes, startZ,
                                    IsCharacter, Rv, IsInteger, myReps, f1, f2, permNonTrivial);
        }
    }
}

