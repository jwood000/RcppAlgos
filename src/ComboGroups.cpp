#include "CleanConvert.h"
#include "NthResult.h"
#include "GmpCombPermUtils.h"
#include "RMatrix.h"
#include <RcppThread.h>

// ******* Overview of the Crucial Part of the Algorithm *******
// -------------------------------------------------------------
// last1 is one plus the upper bound in the previous section, so to obtain the current
// current upper bound, we must first add the size of a section (i.e. grpSize) and su-
// bstract one. We can now compute the length we need to reset v by subtracting idx1. E.g.
//
// Given a portion of v w/ s1 = 9, gSize = 4, idx1 = 9, 6 groups (24 subjects) and base 0:
//       
//              prev sections   bound (index = 8)
//                  /  \        |
//            ............. 8 | 9 12 23 24 | 10 20 21 22 | 11 ... 
//                                 |
//                               idx1 (equal to last1, in this case)
//
// Sort v past idx1:
//                      ... 8 | 9 12 10 11 | 13 14 15 16 | 17... 
//
// Determine the index, idx3, such that v[idx3] > v[idx1]
//
//                      ... 8 | 9 12 10 11 | 13 14 15 16 | 17 ... 
//                                 |          |
//                               idx1       idx3
// 
// Swap idx1 and idx3:
//                      ... 8 | 9 13 10 11 | 12 14 15 16 | 17... 
//
// Move enough indices after idx1 to fill that specific group:
//
//                      ... 8 | 9 13 __ __ | 10 11 12 14 | 15 16 ... 
//
// Identify and move indices that are successively incrementing values of v past idx1:
//
//                      ... 8 | 9 13 14 15 | 10 11 12 16 | 17 ... 
//
// The last two steps are accomplished with std::rotate. This completes the algorithm.

bool nextComboGroup(std::vector<int> &z, int nGrps, 
                    int grpSize, int idx1, int idx2, int last1) {
    
    while (idx2 > idx1 && z[idx2] > z[idx1])
        --idx2;
    
    if ((idx2 + 1) < static_cast<int>(z.size())) {
        if (z[idx2 + 1] > z[idx1])
            std::swap(z[idx1], z[idx2 + 1]);
        
        return true;
    } else {
        const auto zbeg = z.begin();
        
        while (idx1 > 0) {
            const int tipPnt = z[idx2];
            
            while (idx1 > last1 && tipPnt < z[idx1])
                --idx1;
            
            if (tipPnt > z[idx1]) { // **Crucial Part**
                int idx3 = idx1 + 1;
                std::sort(zbeg + idx3, z.end());
                const int xtr = last1 + grpSize - idx3;
                
                while (z[idx3] < z[idx1])
                    ++idx3;
                
                std::swap(z[idx3], z[idx1]);
                std::rotate(zbeg + idx1 + 1, 
                            zbeg + idx3 + 1, zbeg + idx3 + xtr);
                return true;
            } else {
                idx1 -= 2;
                idx2 -= grpSize;
                last1 -= grpSize;
            }
        }
    }
    
    return false;
}

double numGroupCombs(int n, int numGroups, int grpSize) {
    
    double result = 1;
    
    for (double i = n; i > numGroups; --i)
        result *= i;
    
    if (result < std::numeric_limits<double>::max()) {
        double myDiv = 1;
        
        for (double i = 2; i <= grpSize; ++i)
            myDiv *= i;
        
        result /= std::pow(myDiv, numGroups);
        return std::round(result);
    } else {
        return std::numeric_limits<double>::infinity();
    }
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
    const int g = gSize - 1;
    int temp = static_cast<int>(nChooseK(s, g));
    int secLen = total / temp;
    
    std::vector<int> res(n, 0);
    std::vector<int> v(s);
    std::iota(v.begin(), v.end(), 1);
    
    int myMin = 0;
    mpz_t mpzDefault;
    mpz_init(mpzDefault);
    
    for (int j = 0; j < (r - 1); ++j) {
        ind2 = std::floor(ind2 / secLen);
        res[j * gSize] = myMin;
        const std::vector<int> comb = nthComb(s, g, ind2, mpzDefault, v);
        
        for (int k = j * gSize + 1, i = 0; k < ((j + 1) * gSize); ++k, ++i)
            res[k] = v[comb[i]];
        
        v.clear();
        
        for (int i = 1; i <= n; ++i) {
            const auto it = std::find(res.begin(), res.end(), i);
            
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
    const int g = gSize - 1;
    
    mpz_t temp, secLen;
    mpz_init(temp); mpz_init(secLen);
    
    nChooseKGmp(temp, s, g);
    mpz_divexact(secLen, computedRowMpz, temp);
    
    std::vector<int> res(n, 0);
    std::vector<int> v(s);
    std::iota(v.begin(), v.end(), 1);
    
    int myMin = 0;
    constexpr double dblDefault = 0;
    
    for (int j = 0; j < (r - 1); ++j) {
        mpz_tdiv_q(ind2, ind2, secLen);
        res[j * gSize] = myMin;
        const std::vector<int> comb = nthCombGmp(s, g, dblDefault, ind2, v);

        for (int k = j * gSize + 1, i = 0; k < ((j + 1) * gSize); ++k, ++i)
            res[k] = v[comb[i]];

        v.clear();
    
        for (int i = 1; i <= n; ++i) {
            const auto it = std::find(res.begin(), res.end(), i);
    
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


template <typename typeRcpp>
void FinalTouch(typeRcpp &GroupsMat, bool IsArray, int grpSize, int r, int n, 
                int nRows, bool IsNamed, const std::vector<double> &mySample, 
                mpz_t *const myBigSamp, bool IsGmp) {
    
    std::vector<std::string> myColNames(r, "Grp");
    
    for (int j = 0; j < r; ++j)
        myColNames[j] += std::to_string(j + 1);
    
    if (IsArray) {
        Rcpp::CharacterVector rcppCols(r);
        rcppCols = myColNames;
        GroupsMat.attr("dim") = Rcpp::IntegerVector::create(nRows, grpSize, r);
        GroupsMat.attr("dimnames") = Rcpp::List::create(R_NilValue, R_NilValue, rcppCols);
    } else {
        std::vector<std::string> extendedName;
        
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < grpSize; ++j)
                extendedName.push_back(myColNames[i]);
        
        Rcpp::CharacterVector rcppExtended(n);
        rcppExtended = extendedName;
        Rcpp::colnames(GroupsMat) = rcppExtended; 
    }
    
    if (IsNamed) SetSampleNames(IsGmp, nRows, GroupsMat, mySample, myBigSamp);
}

template <typename typeRcpp, typename typeVector>
void SampleWorker(std::size_t n, const typeVector &v, typeRcpp &GroupsMat,
                  int r, int grpSize, std::size_t strtIdx, std::size_t endIdx,
                  bool IsGmp, const std::vector<double> &mySample, mpz_t *const myBigSamp,
                  double computedRows, mpz_t &computedRowMpz) {
    
    if (IsGmp) {
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthComboGroupGmp(n, grpSize, r, myBigSamp[i], computedRowMpz);
            
            for (std::size_t j = 0; j < n; ++j)
                GroupsMat(i, j) = v[z[j]];
        }
    } else {
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthComboGroup(n, grpSize, r, mySample[i], computedRows);
            
            for (std::size_t j = 0; j < n; ++j)
                GroupsMat(i, j) = v[z[j]];
        }
    }
}

template <typename typeRcpp, typename typeVector>
void GroupWorker(std::size_t n, const typeVector &v, typeRcpp &GroupsMat, std::vector<int> z,
                 int r, int grpSize, std::size_t strtIdx, std::size_t endIdx) {
    
    const int idx1 = (r - 1) * grpSize - 1;
    const int idx2 = v.size() - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const std::size_t lastRow = endIdx - 1;
    
    for (std::size_t i = strtIdx; i < lastRow; ++i, nextComboGroup(z, r, grpSize, idx1, idx2, last1))
        for (std::size_t j = 0; j < n; ++j)
            GroupsMat(i, j) = v[z[j]];
    
    // Get last combo group
    for (std::size_t j = 0; j < n; ++j)
        GroupsMat(lastRow, j) = v[z[j]];
}

template <typename typeRcpp, typename typeVector>
void SerialGlue(std::size_t n, const typeVector &v, typeRcpp &GroupsMat, bool IsGmp,
                const std::vector<int> &z, int r, int grpSize, std::size_t nRows,
                bool IsSample, const std::vector<double> &mySample, mpz_t *const myBigSamp,
                double computedRows, mpz_t &computedRowMpz, bool IsArray, bool IsNamed) {
    
    if (IsSample) {
        SampleWorker(n, v, GroupsMat, r, grpSize, 0, nRows,
                     IsGmp, mySample, myBigSamp, computedRows, computedRowMpz);
    } else {
        GroupWorker(n, v, GroupsMat, z, r, grpSize, 0, nRows);
    }
    
    FinalTouch(GroupsMat, IsArray, grpSize, r, n,
               nRows, IsNamed, mySample, myBigSamp, IsGmp);
}

template <typename typeRcpp, typename typeElem>
void ParallelGlue(std::size_t n, const std::vector<typeElem> &v, typeRcpp &GroupsMat, bool IsGmp,
                  const std::vector<int> &z, int r, int grpSize, std::size_t strtIdx, std::size_t endIdx,
                  bool IsSample, const std::vector<double> &mySample, mpz_t *const myBigSamp,
                  double computedRows, mpz_t &computedRowMpz) {
    
    if (IsSample) {
        SampleWorker(n, v, GroupsMat, r, grpSize, strtIdx, endIdx,
                     IsGmp, mySample, myBigSamp, computedRows, computedRowMpz);
    } else {
        GroupWorker(n, v, GroupsMat, z, r, grpSize, strtIdx, endIdx);
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
void GroupsMaster(std::size_t n, const std::vector<typeElem> &v, typeRcpp &GroupsMat, std::vector<int> z,
                  int r, int grpSize, std::size_t nRows, bool Parallel, int nThreads, bool IsGmp,
                  double lower, mpz_t &lowerMpz, double computedRows, mpz_t &computedRowMpz, 
                  bool IsSample, const std::vector<double> &mySample, mpz_t *const myBigSamp,
                  bool IsArray, bool IsNamed) {
    
    if (Parallel) {
        RcppParallel::RMatrix<typeElem> parMat(GroupsMat);
        RcppThread::ThreadPool pool(nThreads);
        std::size_t step = 0, stepSize = nRows / nThreads;
        std::size_t nextStep = stepSize;

        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
            pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), n,
                      std::cref(v), std::ref(parMat), IsGmp, z, r, grpSize, step, nextStep, IsSample, 
                      std::cref(mySample), myBigSamp, computedRows, std::ref(computedRowMpz));
             
             GetStartGrp(IsGmp, z, n, grpSize, r, lower, 
                         computedRows, lowerMpz, computedRowMpz, stepSize);
        }

        pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), n,
                  std::cref(v), std::ref(parMat), IsGmp, z, r, grpSize, step, nRows, IsSample, 
                  std::cref(mySample), myBigSamp, computedRows, std::ref(computedRowMpz));

        pool.join();
        
        FinalTouch(GroupsMat, IsArray, grpSize, r,
                   n, nRows, IsNamed, mySample, myBigSamp, IsGmp);
    } else {
        SerialGlue(n, v, GroupsMat, IsGmp, z, r, grpSize, nRows, IsSample,
                   mySample, myBigSamp, computedRows, computedRowMpz, IsArray, IsNamed);
    }
}

// [[Rcpp::export]]
SEXP ComboGroupsRcpp(SEXP Rv, SEXP RNumGroups, SEXP RRetType, SEXP Rlow, 
                     SEXP Rhigh, bool IsFactor, bool IsCount, SEXP Rparallel,
                     SEXP RNumThreads, int maxThreads, bool IsSample, SEXP RindexVec,
                     SEXP RmySeed, SEXP RNumSamp, Rcpp::Function baseSample, SEXP RNamed) {
    
    int n, numGroups;
    CleanConvert::convertPrimitive(RNumGroups, numGroups, "numGroups");
    bool IsLogical, IsCharacter, IsInteger, IsComplex, IsRaw;
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsNamed = (IsSample) ? CleanConvert::convertLogical(RNamed, "namedSample") : false;
    
    std::vector<double> vNum;
    std::vector<int> vInt;
    Rcpp::CharacterVector rcppChar;
    Rcpp::ComplexVector rcppCplx;
    Rcpp::RawVector rcppRaw;
    
    SetClass(IsCharacter, IsLogical, IsInteger, IsComplex, IsRaw, Rv);
    SetValues(IsCharacter, IsLogical, IsInteger, IsComplex, IsRaw,
              rcppChar, vInt, vNum, rcppCplx, rcppRaw, n, Rv);
    
    if (IsFactor)
        IsCharacter = IsInteger = false;
    
    if (n % numGroups != 0) {
        Rcpp::stop("The length of v (if v is a vector) or v (if v is a"
                       " scalar) must be divisible by numGroups");
    }
    
    const int grpSize = n / numGroups;
    const double computedRows = numGroupCombs(n, numGroups, grpSize);
    bool IsGmp = (computedRows > Significand53);
    
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
    
    if (!IsSample) {
        SetBounds(IsCount, Rlow, Rhigh, IsGmp, bLower, bUpper, 
                  lower, upper, lowerMpz.get(), upperMpz.get());
        CheckBounds(IsGmp, lower, upper, computedRows, lowerMpz[0], upperMpz[0], computedRowMpz);
    }
    
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
    
    int nRows = 0;
    
    if (!IsSample) {
        double userNumRows = 0;
        bool permHolder = false;
        SetNumResults(IsGmp, bLower, bUpper, false, permHolder, upperMpz.get(), 
                      lowerMpz.get(), lower, upper, computedRows, computedRowMpz, nRows, userNumRows);
    }
    
    const std::string retType = Rcpp::as<std::string>(RRetType);
    bool IsArray = false;

    if (retType != "3Darray" && retType != "matrix") {
        Rcpp::stop("retType must be '3Darray' or 'matrix'");
    } else {
        if (retType == "3Darray")
            IsArray = true;
    }
    
    std::size_t sampSize;
    std::vector<double> mySample;
    
    if (IsSample) {
        // sampleLimit defined in GmpCombPermUtils.cpp.. see comments in GmpCombPermUtils.cpp
        IsGmp = computedRows > sampleLimit;
        SetRandomSample(RindexVec, RNumSamp, sampSize, IsGmp, computedRows, mySample, baseSample);
    }
    
    int nThreads = 1;
    const int limit = (IsSample) ? 2 : 20000;
    const int numResults = (IsSample) ? sampSize : nRows;
    SetThreads(Parallel, maxThreads, numResults, IsCharacter, nThreads, RNumThreads, limit);
    
    const std::size_t bigSampSize = (IsSample) ? sampSize : 1;
    auto myVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);
    
    for (std::size_t i = 0; i < bigSampSize; ++i)
        mpz_init(myVec[i]);
        
    if (IsSample) {
        nRows = sampSize;
        SetRandomSampleMpz(RindexVec, RmySeed, sampSize, IsGmp, computedRowMpz, myVec.get());
    }
    
    if (IsCharacter) {
        Rcpp::CharacterMatrix charGroupsMat = Rcpp::no_init_matrix(nRows, n);
        SerialGlue(n, rcppChar, charGroupsMat, IsGmp, startZ, 
                   numGroups, grpSize, nRows, IsSample, mySample, 
                   myVec.get(), computedRows, computedRowMpz, IsArray, IsNamed);
        return charGroupsMat;
    } else if (IsComplex) {
        Rcpp::ComplexMatrix cplxGroupsMat = Rcpp::no_init_matrix(nRows, n);
        SerialGlue(n, rcppCplx, cplxGroupsMat, IsGmp, startZ, 
                   numGroups, grpSize, nRows, IsSample, mySample, 
                   myVec.get(), computedRows, computedRowMpz, IsArray, IsNamed);
        return cplxGroupsMat;
    } else if (IsRaw) {
        Rcpp::RawMatrix rawGroupsMat = Rcpp::no_init_matrix(nRows, n);
        SerialGlue(n, rcppRaw, rawGroupsMat, IsGmp, startZ, 
                   numGroups, grpSize, nRows, IsSample, mySample, 
                   myVec.get(), computedRows, computedRowMpz, IsArray, IsNamed);
        return rawGroupsMat;
    } else if (IsLogical) {
        Rcpp::LogicalMatrix boolGroupsMat(nRows, n);
        GroupsMaster(n, vInt, boolGroupsMat, startZ, numGroups, grpSize, nRows, Parallel,
                     nThreads, IsGmp, lower, lowerMpz[0], computedRows, computedRowMpz,
                     IsSample, mySample, myVec.get(), IsArray, IsNamed);
        return boolGroupsMat;
    } else if (IsInteger || IsFactor) {
        Rcpp::IntegerMatrix intGroupsMat = Rcpp::no_init_matrix(nRows, n);
        GroupsMaster(n, vInt, intGroupsMat, startZ, numGroups, grpSize, nRows, Parallel,
                     nThreads, IsGmp, lower, lowerMpz[0], computedRows, computedRowMpz,
                     IsSample, mySample, myVec.get(), IsArray, IsNamed);

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
                     nThreads, IsGmp, lower, lowerMpz[0], computedRows, computedRowMpz,
                     IsSample, mySample, myVec.get(), IsArray, IsNamed);
        return numGroupsMat;
    }
}
