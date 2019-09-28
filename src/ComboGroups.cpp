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
void GroupWorker(std::size_t n, typeVector v, typeRcpp &GroupsMat, std::vector<int> z,
                 int r, int grpSize, std::size_t strtIdx, std::size_t endIdx) {
    
    const int idx1 = (r - 1) * grpSize - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const std::size_t lastRow = endIdx - 1;
    
    for (std::size_t i = strtIdx; i < lastRow; ++i, nextComboGroup(z, r, grpSize, idx1, last1))
        for (std::size_t j = 0; j < n; ++j)
            GroupsMat(i, j) = v[z[j]];
        
    // Get last combo group
    for (std::size_t j = 0; j < n; ++j)
        GroupsMat(lastRow, j) = v[z[j]];
}

template <typename typeRcpp, typename typeElem>
void ParallelGlue(std::size_t n, std::vector<typeElem> v, typeRcpp &GroupsMat,
                  std::vector<int> z, int r, int grpSize, std::size_t strtIdx, std::size_t endIdx) {
    GroupWorker(n, v, GroupsMat, z, r, grpSize, strtIdx, endIdx);
}

template <typename typeRcpp>
void FinalTouch(typeRcpp &GroupMat, bool IsArray,
                int grpSize, int r, int n, int nRows) {
    
    std::vector<std::string> myColNames(r, "Grp");
    
    for (int j = 0; j < r; ++j)
        myColNames[j] += std::to_string(j + 1);
    
    if (IsArray) {
        Rcpp::CharacterVector rcppCols(r);
        rcppCols = myColNames;
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

void GetStartGrp(bool IsGmp, std::vector<int> &z, int n, int grpSize, 
                 int r, double &lower, double computedRows, mpz_t lowerMpz,
                 mpz_t computedRowMpz, std::size_t stepSize) {
    if (IsGmp) {
        mpz_add_ui(lowerMpz, lowerMpz, stepSize);
        z = nthComboGroupGmp(n, grpSize, r, lowerMpz, computedRowMpz);
    } else {
        lower += stepSize;
        z = nthComboGroup(n, grpSize, r, lower, computedRows);
    }
}

template <typename typeRcpp, typename typeElem>
void GroupsMaster(std::size_t n, std::vector<typeElem> v, typeRcpp &GroupsMat, std::vector<int> z,
                  int r, int grpSize, std::size_t nRows, bool Parallel, int nThreads, bool IsGmp,
                  double lower, mpz_t &lowerMpz, double computedRows, mpz_t &computedRowMpz) {
    
    if (Parallel) {
        RcppParallel::RMatrix<typeElem> parMat(GroupsMat);
        RcppThread::ThreadPool pool(nThreads);
        std::size_t step = 0, stepSize = nRows / nThreads;
        std::size_t nextStep = stepSize;

        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
            pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), n,
                      v, std::ref(parMat), z, r, grpSize, step, nextStep);
             
             GetStartGrp(IsGmp, z, n, grpSize, r, lower, 
                         computedRows, lowerMpz, computedRowMpz, stepSize);
        }

        pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), n,
                  v, std::ref(parMat), z, r, grpSize, step, nRows);

        pool.join();
    } else {
        GroupWorker(n, v, GroupsMat, z, r, grpSize, 0, nRows);
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
    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);
    
    mpz_init(lowerMpz[0]); mpz_init(upperMpz[0]);
    SetBounds(IsCount, Rlow, Rhigh, IsGmp, bLower, bUpper, 
              lower, upper, lowerMpz.get(), upperMpz.get());
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
    SetNumResults(IsGmp, bLower, bUpper, false, permHolder, upperMpz.get(), 
                  lowerMpz.get(), lower, upper, computedRows, computedRowMpz, nRows, userNumRows);

    const std::string retType = Rcpp::as<std::string>(RRetType);
    bool isArray = false;

    if (retType != "3Darray" && retType != "matrix") {
        Rcpp::stop("retType must be '3Darray' or 'matrix'");
    } else {
        if (retType == "3Darray")
            isArray = true;
    }

    int nThreads = 1;
    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows, IsCharacter, nThreads, RNumThreads, limit);

    if (IsCharacter) {
        Rcpp::CharacterMatrix charGroupsMat = Rcpp::no_init_matrix(nRows, n);
        GroupWorker(n, rcppChar, charGroupsMat, startZ, numGroups, grpSize, 0u, nRows);
        FinalTouch(charGroupsMat, isArray, grpSize, numGroups, n, nRows);
        return charGroupsMat;
    } else if (IsLogical) {
        Rcpp::LogicalMatrix boolGroupsMat(nRows, n);
        GroupsMaster(n, vInt, boolGroupsMat, startZ, numGroups, grpSize, nRows, Parallel, 
                     nThreads, IsGmp, lower, lowerMpz[0], computedRows, computedRowMpz);
        FinalTouch(boolGroupsMat, isArray, grpSize, numGroups, n, nRows);
        return boolGroupsMat;
    } else if (IsInteger || IsFactor) {
        Rcpp::IntegerMatrix intGroupsMat = Rcpp::no_init_matrix(nRows, n);
        GroupsMaster(n, vInt, intGroupsMat, startZ, numGroups, grpSize, nRows, Parallel, 
                     nThreads, IsGmp, lower, lowerMpz[0], computedRows, computedRowMpz);
        FinalTouch(intGroupsMat, isArray, grpSize, numGroups, n, nRows);

        if (IsFactor) {
            Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
            Rcpp::CharacterVector myClass = testFactor.attr("class");
            Rcpp::CharacterVector myLevels = testFactor.attr("levels");
            intGroupsMat.attr("class") = myClass;
            intGroupsMat.attr("levels") = myLevels;
        }

        return intGroupsMat;
    } else {
        Rcpp::NumericMatrix numGroupsMat = Rcpp::no_init_matrix(nRows, n);
        GroupsMaster(n, vNum, numGroupsMat, startZ, numGroups, grpSize, nRows, Parallel, 
                     nThreads, IsGmp, lower, lowerMpz[0], computedRows, computedRowMpz);
        FinalTouch(numGroupsMat, isArray, grpSize, numGroups, n, nRows);
        return numGroupsMat;
    }
}
