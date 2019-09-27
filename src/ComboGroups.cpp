#include "CleanConvert.h"
#include "NthResult.h"
#include "CountGmp.h"
#include "RMatrix.h"
#include <RcppThread.h>

// ******* Overview of the Crucial Part of the Algorithm *******
// -------------------------------------------------------------
// last1 is one plus the upper bound in the previous section, so
// to obtain the current upper bound, we must first add the size
// of a section (i.e. gSize) and substract one. We can now comp-
// ute the length we need to reset v by subtracting idx1. E.g.
//
// Given a portion of v w/ s1 = 9, gSize = 4, idx1 = 9,
//       6 groups (24 total subjects) and base 0:
//
//  prev sections   bound (index = 8)
//        /  \        |
//  ............. 8 | 9 12 23 24 | 10 20 21 22 | 11 ... 
//                       |
//                     idx1 (equal to last1, in this case)
//
// Sort v past idx1:
//            ... 8 | 9 12 10 11 | 13 14 15 16 | 17... 
//
// Determine the index, idx3, such that v[idx3] > v[idx1]
//
//            ... 8 | 9 12 10 11 | 13 14 15 16 | 17 ... 
//                       |          |
//                      idx1       idx3
// 
// Swap idx1 and idx3:
//            ... 8 | 9 13 10 11 | 12 14 15 16 | 17... 
//
// Move enough indices after idx1 to fill that specific group:
//
//            ... 8 | 9 13 __ __ | 10 11 12 14 | 15 16 ... 
//
// Identify and move indices that are successively incrementing
// values of v past idx1:
//
//            ... 8 | 9 13 14 15 | 10 11 12 16 | 17 ... 
//
// The last two steps are accomplished with std::rotate. This
// completes the algorithm.

bool nextComboGroup(std::vector<int> &v, int nGrps, 
                    int grpSize, int idx1, int last1) {
    
    int idx2 = v.size() - 1;
    int last2 = idx1;
    const int lastIdx = idx2;
    
    while (idx2 > last2 && v[idx2] > v[idx1])
        --idx2;
    
    if (idx2 < lastIdx) {
        if (v[idx2 + 1] > v[idx1])
            std::swap(v[idx1], v[idx2 + 1]);
        
        return true;
    } else {
        while (idx1 > 0) {
            while (idx1 > last1 && v[idx2] < v[idx1])
                --idx1;
            
            if (v[idx2] > v[idx1]) { // **Crucial Part**
                int idx3 = idx1 + 1;
                const int lenDiff = last1 + grpSize - 1 - idx1;
                std::sort(v.begin() + idx3, v.end());
                
                while (v[idx3] < v[idx1])
                    ++idx3;
                
                std::swap(v[idx3], v[idx1]);
                std::rotate(v.begin() + idx1 + 1, v.begin() 
                                + idx3 + 1, v.begin() + idx3 + lenDiff);
                
                return true;
            } else {
                idx1 -= 2;
                last1 -= grpSize;
                idx2 -= grpSize;
            }
        }
    }
    
    return false;
}

double numGroupCombs(int n, int numGroups, int grpSize) {
    
    double result = 1;
    
    for (int i = n; i > numGroups; --i)
        result *= i;
    
    double myDiv = 1;
    
    for (int i = 2; i <= grpSize; ++i)
        myDiv *= i;
    
    result /= std::pow(myDiv, numGroups);
    return std::round(result);
}

void numGroupCombsGmp(mpz_t result, int n, 
                      int numGroups, int grpSize) {

    for (int i = n; i > numGroups; --i)
        mpz_mul_ui(result, result, i);

    mpz_t myDiv;
    mpz_init(myDiv);
    mpz_set_ui(myDiv, 1);

    for (int i = 2; i <= grpSize; ++i)
        mpz_mul_ui(myDiv, myDiv, i);
    
    mpz_pow_ui(myDiv, myDiv, numGroups);
    mpz_divexact(result, result, myDiv);
    mpz_clear(myDiv);
}

std::vector<int> nthComboGroup(int n, int gSize, int r,
                               double myIndex, double total) {
    
    double ind1 = myIndex, ind2 = myIndex;
    int s = n - 1;
    int g = gSize - 1;
    int temp = static_cast<int>(nChooseK(s, g));
    int secLen = total / temp;
    
    std::vector<int> res(n, 0);
    std::vector<int> v(s);
    std::iota(v.begin(), v.end(), 1);
    std::vector<int>::iterator it;
    
    int myMin = 0;
    mpz_t mpzDefault;
    mpz_init(mpzDefault);
    
    for (int j = 0; j < (r - 1); ++j) {
        ind2 = std::floor(ind2 / secLen);
        res[j * gSize] = myMin;
        std::vector<int> comb = nthComb(s, g, ind2, mpzDefault, v);
        
        for (int k = j * gSize + 1, i = 0; k < ((j + 1) * gSize); ++k, ++i)
            res[k] = v[comb[i]];
        
        v.clear();
        
        for (int i = 1; i <= n; ++i) {
            it = std::find(res.begin(), res.end(), i);
            
            if (it == res.end())
                v.push_back(i);
        }
        
        myMin = v.front();
        v.erase(v.begin());
        ind1 -= ind2 * secLen;
        ind2 = ind1;
        s -= gSize;
        temp = static_cast<int>(nChooseK(s, g));
        secLen /= temp;
    }
    
    res[(r - 1) * gSize] = myMin;
    
    for (int k = (r - 1) * gSize + 1, i = 0; k < (r * gSize); ++k, ++i)
        res[k] = v[i];
    
    return res;
}

std::vector<int> nthComboGroupGmp(int n, int gSize, int r,
                                  mpz_t lowerMpz, mpz_t computedRowMpz) {
    mpz_t ind1, ind2;
    mpz_init(ind1); mpz_init(ind2);
    mpz_set(ind1, lowerMpz); mpz_set(ind2, lowerMpz);
    
    int s = n - 1;
    int g = gSize - 1;
    
    mpz_t temp, secLen;
    mpz_init(temp); mpz_init(secLen);
    
    nChooseKGmp(temp, s, g);
    mpz_divexact(secLen, computedRowMpz, temp);
    
    std::vector<int> res(n, 0);
    std::vector<int> v(s);
    std::iota(v.begin(), v.end(), 1);
    std::vector<int>::iterator it;
    
    int myMin = 0;
    const double dblDefault = 0;
    
    for (int j = 0; j < (r - 1); ++j) {
        mpz_tdiv_q(ind2, ind2, secLen);
        res[j * gSize] = myMin;
        std::vector<int> comb = nthCombGmp(s, g, dblDefault, ind2, v);

        for (int k = j * gSize + 1, i = 0; k < ((j + 1) * gSize); ++k, ++i)
            res[k] = v[comb[i]];

        v.clear();
    
        for (int i = 1; i <= n; ++i) {
            it = std::find(res.begin(), res.end(), i);
    
            if (it == res.end())
                v.push_back(i);
        }
    
        myMin = v.front();
        v.erase(v.begin());
        mpz_mul(temp, ind2, secLen);
        mpz_sub(ind1, ind1, temp);
        mpz_set(ind2, ind1);
    
        s -= gSize;
        nChooseKGmp(temp, s, g);
        mpz_divexact(secLen, secLen, temp);
    }

    res[(r - 1) * gSize] = myMin;

    for (int k = (r - 1) * gSize + 1, i = 0; k < (r * gSize); ++k, ++i)
        res[k] = v[i];
    
    return res;
}

template <typename typeRcpp, typename typeVector>
void GroupMaster(std::size_t n, typeVector v, typeRcpp GroupMat, Rcpp::List myList,
                 std::vector<int> z, std::size_t nRows, int r, int grpSize, bool isArray,
                 bool isList, SEXP Rv, bool IsFactor) {
    
    const int idx1 = (r - 1) * grpSize - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const int lastRow = nRows - 1;
    std::vector<std::string> myColNames(r, "Grp");
    
    for (int j = 0; j < r; ++j)
        myColNames[j] += std::to_string(j + 1);
    
    Rcpp::CharacterVector rcppCols(myColNames.size());
    rcppCols = myColNames;
    
    if (isList) {
        if (IsFactor) {
            Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
            Rcpp::CharacterVector myClass = testFactor.attr("class");
            Rcpp::CharacterVector myLevels = testFactor.attr("levels");
            
            for (std::size_t i = 0; i < lastRow; ++i, nextComboGroup(z, r, grpSize, idx1, last1)) {
                typeRcpp tempMat(grpSize, r);
                for (int j = 0, ind = 0; j < r; ++j)
                    for (std::size_t k = 0; k < grpSize; ++k, ++ind)
                        tempMat(k, j) = v[z[ind]];
                
                tempMat.attr("class") = myClass;
                tempMat.attr("levels") = myLevels;
                Rcpp::colnames(tempMat) = rcppCols;
                myList[i] = tempMat;
            }
            
            // Get last combo group
            for (int j = 0, ind = 0; j < r; ++j)
                for (std::size_t k = 0; k < grpSize; ++k, ++ind)
                    GroupMat(k, j) = v[z[ind]];
            
            GroupMat.attr("class") = myClass;
            GroupMat.attr("levels") = myLevels;
            Rcpp::colnames(GroupMat) = rcppCols;
            myList[lastRow] = GroupMat;
        } else {
            for (std::size_t i = 0; i < lastRow; ++i, nextComboGroup(z, r, grpSize, idx1, last1)) {
                typeRcpp tempMat(grpSize, r);
                
                for (int j = 0, ind = 0; j < r; ++j)
                    for (std::size_t k = 0; k < grpSize; ++k, ++ind)
                        tempMat(k, j) = v[z[ind]];
                
                Rcpp::colnames(tempMat) = rcppCols;
                myList[i] = tempMat;
            }
            
            // Get last combo group
            for (int j = 0, ind = 0; j < r; ++j)
                for (std::size_t k = 0; k < grpSize; ++k, ++ind)
                    GroupMat(k, j) = v[z[ind]];
            
            Rcpp::colnames(GroupMat) = rcppCols;
            myList[lastRow] = GroupMat;
        }
    } else {
        for (std::size_t i = 0; i < lastRow; ++i, nextComboGroup(z, r, grpSize, idx1, last1))
            for (std::size_t j = 0; j < n; ++j)
                GroupMat(i, j) = v[z[j]];
            
        // Get last combo group
        for (std::size_t j = 0; j < n; ++j)
            GroupMat(lastRow, j) = v[z[j]];
        
        if (isArray) {
            GroupMat.attr("dim") = Rcpp::IntegerVector::create(nRows, grpSize, r);
            GroupMat.attr("dimnames") = Rcpp::List::create(R_NilValue, R_NilValue, rcppCols);
        } else {
            std::vector<std::string> extendedName;
            
            for (int i = 0; i < r; ++i)
                for (int j = 0; j < grpSize; ++j)
                    extendedName.push_back(myColNames[i]);
            
            Rcpp::CharacterVector rcppExtended(n);
            rcppExtended = extendedName;
            Rcpp::colnames(GroupMat) = rcppExtended;
        }
    }
}

// [[Rcpp::export]]
SEXP ComboGroupsRcpp(SEXP Rv, SEXP RNumGroups, SEXP RRetType, SEXP Rlow, 
                    SEXP Rhigh, bool IsFactor, bool IsCount, SEXP Rparallel,
                    SEXP RNumThreads, int maxThreads) {
    
    int n, numGroups;
    CleanConvert::convertPrimitive(RNumGroups, numGroups, "numGroups");
    bool IsLogical, IsCharacter, IsInteger;
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    
    std::vector<double> vNum;
    std::vector<int> vInt;
    Rcpp::CharacterVector rcppChar;
    
    switch(TYPEOF(Rv)) {
        case LGLSXP: {
            IsLogical = true;
            IsInteger = IsCharacter = false;
            break;
        }
        case INTSXP: {
            IsInteger = true;
            IsLogical = IsCharacter = false;
            break;
        }
        case REALSXP: {
            IsLogical = IsInteger = IsCharacter = false;
            break;
        }
        case STRSXP: {
            IsCharacter = true;
            Parallel = IsLogical = IsInteger = false;
            break;
        }
        default: {
            Rcpp::stop("Only integers, numerical, character, and factor classes are supported for v");   
        }
    }
    
    SetValues(IsCharacter, IsLogical, IsInteger, rcppChar, vInt, vNum, n, Rv);
    
    if (IsFactor)
        IsCharacter = IsInteger = false;
    
    if (n % numGroups != 0) {
        Rcpp::stop("The length of v (if v is a vector) or v (if v is a"
                       " scalar) must be divisible by numGroups");
    }
    
    const int grpSize = n / numGroups;
    const double computedRows = numGroupCombs(n, numGroups, grpSize);
    const bool IsGmp = (computedRows > Significand53);
    
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);
    
    if (IsGmp) {
        mpz_set_ui(computedRowMpz, 1);
        numGroupCombsGmp(computedRowMpz, n, numGroups, grpSize);
    }
    
    double lower = 0, upper = 0;
    bool bLower = false, bUpper = false;
    mpz_t lowerMpz[1], upperMpz[1];
    mpz_init(lowerMpz[0]); mpz_init(upperMpz[0]);
    mpz_set_ui(lowerMpz[0], 0); mpz_set_ui(upperMpz[0], 0);
    
    if (!IsCount) {
        if (!Rf_isNull(Rlow)) {
            bLower = true;
            
            if (IsGmp) {
                createMPZArray(Rlow, lowerMpz, 1, "lower");
                mpz_sub_ui(lowerMpz[0], lowerMpz[0], 1);
            } else {                                    // numOnly = false
                CleanConvert::convertPrimitive(Rlow, lower, "lower", false);
                --lower;
            }
        }
        
        if (!Rf_isNull(Rhigh)) {
            bUpper = true;
            
            if (IsGmp) {
                createMPZArray(Rhigh, upperMpz, 1, "upper");
            } else {                                     // numOnly = false
                CleanConvert::convertPrimitive(Rhigh, upper, "upper", false);
            }
        }
    }
    
    CheckBounds(IsGmp, lower, upper, computedRows, lowerMpz[0], upperMpz[0], computedRowMpz);
    
    if (IsCount)
        return GetCount(IsGmp, computedRowMpz, computedRows);
    
    std::vector<int> startZ;
    double dblLower = lower;
    if (!IsGmp) mpz_set_d(lowerMpz[0], dblLower);
    
    if (bLower && mpz_cmp_ui(lowerMpz[0], 0) > 0) {
        if (IsGmp)
            startZ = nthComboGroupGmp(n, grpSize, numGroups, lowerMpz[0], computedRowMpz);
        else
            startZ = nthComboGroup(n, grpSize, numGroups, lower, computedRows);
    } else {
        startZ.resize(n);
        std::iota(startZ.begin(), startZ.end(), 0);
    }
    
    double userNumRows = 0;
    int nRows = 0;
    bool permHolder = false;
    SetNumResults(IsGmp, bLower, bUpper, false, permHolder, upperMpz, lowerMpz,
                  lower, upper, computedRows, computedRowMpz, nRows, userNumRows);

    const std::string retType = Rcpp::as<std::string>(RRetType);
    bool isArray = false;
    bool isList = false;

    if (retType != "3Darray" && retType != "list" && retType != "matrix") {
        Rcpp::stop("retType must be '3Darray', 'list' or matrix'");
    } else {
        if (retType == "list")
            isList = true;
        else if (retType == "3Darray")
            isArray = true;
    }

    int nThreads = 1;
    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows, IsCharacter, nThreads, RNumThreads, limit);

    if (isList) {
        Rcpp::List myList(nRows);
        if (IsCharacter) {
            Rcpp::CharacterMatrix charGrpMat(grpSize, numGroups);
            GroupMaster(n, rcppChar, charGrpMat, myList, startZ, nRows,
                        numGroups, grpSize, isArray, isList, Rv, IsFactor);
        } else if (IsInteger || IsFactor) {
            Rcpp::IntegerMatrix intGrpMat(grpSize, numGroups);
            GroupMaster(n, vInt, intGrpMat, myList, startZ, nRows,
                        numGroups, grpSize, isArray, isList, Rv, IsFactor);
        } else {
            Rcpp::NumericMatrix numGrpMat(grpSize, numGroups);
            GroupMaster(n, vNum, numGrpMat, myList, startZ, nRows,
                        numGroups, grpSize, isArray, isList, Rv, IsFactor);
        }

        return myList;
    } else {
        Rcpp::List TrivialList;
        if (IsCharacter) {
            Rcpp::CharacterMatrix charGroupMat = Rcpp::no_init_matrix(nRows, n);
            GroupMaster(n, rcppChar, charGroupMat, TrivialList, startZ,
                        nRows, numGroups, grpSize, isArray, isList, Rv, IsFactor);
            return charGroupMat;
        } else if (IsInteger || IsFactor) {
            Rcpp::IntegerMatrix intGroupMat = Rcpp::no_init_matrix(nRows, n);
            GroupMaster(n, vInt, intGroupMat, TrivialList, startZ,
                        nRows, numGroups, grpSize, isArray, isList, Rv, IsFactor);

            if (IsFactor) {
                Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
                Rcpp::CharacterVector myClass = testFactor.attr("class");
                Rcpp::CharacterVector myLevels = testFactor.attr("levels");
                intGroupMat.attr("class") = myClass;
                intGroupMat.attr("levels") = myLevels;
            }

            return intGroupMat;
        } else {
            Rcpp::NumericMatrix numGroupMat = Rcpp::no_init_matrix(nRows, n);
            GroupMaster(n, vNum, numGroupMat, TrivialList, startZ,
                        nRows, numGroups, grpSize, isArray, isList, Rv, IsFactor);
            return numGroupMat;
        }
    }
}
