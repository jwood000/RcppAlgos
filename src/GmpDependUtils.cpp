#include "GmpDependUtils.h"
#include "ConstraintsUtils.h"
#include "PartitionsCounts.h"
#include "CheckStdRet.h"

static gmp_randstate_t seed_state;
static int seed_init = 0;

constexpr std::size_t intSize = sizeof(int);

SEXP GetCount(bool IsGmp, mpz_t computedRowsMpz, double computedRows) {
    
    if (IsGmp) {
        constexpr std::size_t numb = 8 * intSize;
        const std::size_t sizeNum = intSize * 
                (2 + (mpz_sizeinbase(computedRowsMpz, 2) + numb - 1) / numb);
        const std::size_t size = intSize + sizeNum;
        
        Rcpp::RawVector ansPos = Rcpp::no_init_vector(size);
        char* rPos = (char*) RAW(ansPos);
        ((int*) rPos)[0] = 1; // first int is vector-size-header

        // current position in rPos[] (starting after vector-size-header)
        myRaw(&rPos[intSize], computedRowsMpz, sizeNum);
        ansPos.attr("class") = Rcpp::CharacterVector::create("bigz");
        return(ansPos);
    } else {
        if (computedRows > std::numeric_limits<int>::max()) {
            return Rcpp::wrap(computedRows);
        } else {
            return Rcpp::wrap(static_cast<int>(computedRows));
        }
    }
}

void SetBounds(const SEXP &Rlow, const SEXP &Rhigh, bool IsGmp, bool &bLower, bool &bUpper,
               double &lower, double &upper, mpz_t *const lowerMpz, 
               mpz_t *const upperMpz, mpz_t computedRowsMpz, double computedRows) {
    
    if (!Rf_isNull(Rlow)) {
        if (IsGmp) {
            createMPZArray(Rlow, lowerMpz, 1, "lower");
            bLower = mpz_cmp_si(lowerMpz[0], 1) > 0;
            lower = (bLower) ? 1 : 0;
            
            if (mpz_cmp(lowerMpz[0], computedRowsMpz) > 0)
                Rcpp::stop("bounds cannot exceed the maximum number of possible results");
            
            mpz_sub_ui(lowerMpz[0], lowerMpz[0], 1);
        } else {                                    // numOnly = false
            CleanConvert::convertPrimitive(Rlow, lower, "lower", false);
            bLower = lower > 1;
            
            if (lower > computedRows)
                Rcpp::stop("bounds cannot exceed the maximum number of possible results");
            
            --lower;
        }
    }
    
    if (!Rf_isNull(Rhigh)) {
        bUpper = true;
        
        if (IsGmp) {
            createMPZArray(Rhigh, upperMpz, 1, "upper");
            
            if (mpz_cmp(upperMpz[0], computedRowsMpz) > 0)
                Rcpp::stop("bounds cannot exceed the maximum number of possible results");
            
        } else {                                     // numOnly = false
            CleanConvert::convertPrimitive(Rhigh, upper, "upper", false);
            
            if (upper > computedRows)
                Rcpp::stop("bounds cannot exceed the maximum number of possible results");
        }
    }
}

void SetNumResults(bool IsGmp, bool bLower, bool bUpper, bool IsGenCnstrd, mpz_t *const upperMpz,
                   mpz_t *const lowerMpz, double lower, double upper, double computedRows, 
                   mpz_t &computedRowsMpz, int &nRows, double &userNumRows) {
    
    if (IsGmp) {
        mpz_t testBound;
        mpz_init(testBound);
        
        if (bLower && bUpper) {
            mpz_sub(testBound, upperMpz[0], lowerMpz[0]);
            mpz_t absTestBound;
            mpz_init(absTestBound);
            mpz_abs(absTestBound, testBound);
            
            if (mpz_cmp_ui(absTestBound, std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(testBound);
            mpz_clear(absTestBound);
        } else if (bUpper) {
            if (mpz_cmp_d(upperMpz[0], std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(upperMpz[0]);
        } else if (bLower) {
            mpz_sub(testBound, computedRowsMpz, lowerMpz[0]);
            mpz_abs(testBound, testBound);
            
            if (mpz_cmp_d(testBound, std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(testBound);
        }
        
        mpz_clear(testBound);
    } else {
        if (bLower && bUpper)
            userNumRows = upper - lower;
        else if (bUpper)
            userNumRows = upper;
        else if (bLower)
            userNumRows = computedRows - lower;
    }
    
    if (userNumRows == 0) {
        if (bLower && bUpper) {
            // Since lower is decremented and upper isn't, this implies that upper - lower = 0
            // which means that lower is one larger than upper as put in by the user
            
            Rcpp::stop("The number of rows must be positive. Either the lowerBound "
                           "exceeds the maximum number of possible results or the "
                           "lowerBound is greater than the upperBound.");
        } else {
            if (computedRows > std::numeric_limits<int>::max() && !IsGenCnstrd)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = computedRows;
            
            if (!IsGenCnstrd)
                nRows = static_cast<int>(computedRows);
        }
    } else if (userNumRows < 0) {
        Rcpp::stop("The number of rows must be positive. Either the lowerBound "
                       "exceeds the maximum number of possible results or the "
                       "lowerBound is greater than the upperBound.");
    } else if (userNumRows > std::numeric_limits<int>::max()) {
        Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
    } else {
        nRows = static_cast<int>(userNumRows);
    }
}

void SetRandomSampleMpz(const SEXP &RindexVec, const SEXP &RmySeed, std::size_t sampSize,
                        bool IsGmp, mpz_t &computedRowsMpz, mpz_t *const myVec) {
    
    if (IsGmp) {
        if (!Rf_isNull(RindexVec)) {
            createMPZArray(RindexVec, myVec, sampSize, "sampleVec");
            
            // get zero base
            for (std::size_t i = 0; i < sampSize; ++i)
                mpz_sub_ui(myVec[i], myVec[i], 1);
            
        } else {
            // The following code is very similar to the source
            // code of gmp::urand.bigz. The main difference is
            // the use of mpz_urandomm instead of mpz_urandomb
            if (seed_init == 0)
                gmp_randinit_default(seed_state);
            
            seed_init = 1;
            
            if (!Rf_isNull(RmySeed)) {
                mpz_t mpzSeed[1];
                mpz_init(mpzSeed[0]);
                createMPZArray(RmySeed, mpzSeed, 1, "seed");
                gmp_randseed(seed_state, mpzSeed[0]);
                mpz_clear(mpzSeed[0]);
            }
            
            // random number is between 0 and gmpRows[0] - 1
            // so we need to add 1 to each element
            for (std::size_t i = 0; i < sampSize; ++i) {
                mpz_init(myVec[i]);
                mpz_urandomm(myVec[i], seed_state, computedRowsMpz);
            }
        }
        
        mpz_t maxGmp;
        mpz_init(maxGmp);
        mpz_set(maxGmp, myVec[0]);
        
        for (std::size_t i = 1; i < sampSize; ++i)
            if (mpz_cmp(myVec[i], maxGmp) > 0)
                mpz_set(maxGmp, myVec[i]);
            
        if (mpz_cmp(maxGmp, computedRowsMpz) >= 0) {
            Rcpp::stop("One or more of the requested values in sampleVec "
                           "exceeds the maximum number of possible results");
        }
    }
}

void SetStartZ(int n, int m, double &lower, int stepSize, mpz_t &lowerMpz, bool IsRep,
               bool IsComb, bool IsMult, bool IsGmp, const std::vector<int> &myReps,
               const std::vector<int> &freqs, std::vector<int> &z, 
               const nthResultPtr nthResFun) {
    
    if (lower > 0 || stepSize > 0) {
        if (stepSize) {
            if (IsGmp) {
                mpz_add_ui(lowerMpz, lowerMpz, stepSize);
            } else {
                lower += stepSize;
            }
        }
        
        z = nthResFun(n, m, lower, lowerMpz, myReps);
        
        if (!IsComb) {
            if (IsMult) {
                std::vector<int> f(n, 0);
                
                for (int i = 0; i < m; ++i)
                    ++f[z[i]];
                
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < (myReps[i] - f[i]); ++j)
                        z.push_back(i); 
                
            } else if (!IsRep) {
                if (m < n)
                    for (int i = 0; i < n; ++i)
                        if (std::find(z.begin(), z.end(), i) == z.end())
                            z.push_back(i);
            }
        }
    } else {
        if (IsComb) {
            if (IsMult)
                z.assign(freqs.cbegin(), freqs.cbegin() + m);
            else if (IsRep)
                std::fill(z.begin(), z.end(), 0);
            else
                std::iota(z.begin(), z.end(), 0);
        } else {
            if (IsMult) {
                z = freqs;
            } else if (IsRep) {
                std::fill(z.begin(), z.end(), 0);
            } else {
                z.resize(n);
                std::iota(z.begin(), z.end(), 0);
            }
        }
    }
}

// [[Rcpp::export]]
Rcpp::List GetClassVals(bool IsStdRet, SEXP Rv, SEXP Rm, SEXP RisRep,
                        SEXP RFreqs, bool IsComb, bool IsFactor, SEXP stdFun) {
    
    int n, m = 0, lenFreqs = 0;
    VecType myType = VecType::Integer;
    bool IsMult = false;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqs;
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    
    SetType(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    SetFreqsAndM(RFreqs, IsMult, myReps, IsRep, lenFreqs, freqs, Rm, n, m);
    
    const SEXP sexpVec = CopyRv(Rv, vInt, vNum, myType, IsFactor);
    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep, n,
                                                m, Rm, lenFreqs, freqs, myReps);
    
    const bool IsGmp = (computedRows > sampleLimit);
    
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);
    
    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult, 
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }
    
    Rcpp::List freqInfo = Rcpp::List::create(Rcpp::wrap(freqs),
                                             Rcpp::wrap(myReps));
    
    const SEXP sexpNumRows = GetCount(IsGmp, computedRowsMpz, computedRows);
    
    // Needed to determine if nextFullPerm or nextPerm will be called
    const bool IsFullPerm = (IsComb || IsRep) ? false :
                            (m == n || m == static_cast<int>(freqs.size()));
    
    Rcpp::LogicalVector bVec = Rcpp::LogicalVector::create(IsFactor, IsComb, IsMult,
                                                           IsRep, IsGmp, IsFullPerm);
    
    const bool applyFun = !Rf_isNull(stdFun) && !IsFactor;
    
    if (applyFun) {
        if (!Rf_isFunction(stdFun))
            Rcpp::stop("FUN must be a function!");
    }
    
    return Rcpp::List::create(Rcpp::Named("Rv") = sexpVec,
                              Rcpp::Named("Rm") = m,
                              Rcpp::Named("nRows") = sexpNumRows,
                              Rcpp::Named("bVec") = bVec,
                              Rcpp::Named("freqInfo") = freqInfo,
                              Rcpp::Named("applyFun") = applyFun);
}
