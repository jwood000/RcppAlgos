#include "CombPermUtils.h"
#include "CleanConvert.h"
#include "Cpp14MakeUnique.h"
#include <gmp.h>

static gmp_randstate_t seed_state;
static int seed_init = 0;

SEXP GetCount(bool IsGmp, mpz_t computedRowMpz, double computedRows) {
    
    if (IsGmp) {
        constexpr std::size_t numb = 8 * sizeof(int);
        const std::size_t sizeNum = sizeof(int) * 
                (2 + (mpz_sizeinbase(computedRowMpz, 2) + numb - 1) / numb);
        const std::size_t size = sizeof(int) + sizeNum;
        
        Rcpp::RawVector ansPos = Rcpp::no_init_vector(size);
        char* rPos = (char*)(RAW(ansPos));
        ((int*)(rPos))[0] = 1; // first int is vector-size-header

        // current position in rPos[] (starting after vector-size-header)
        std::size_t posPos = sizeof(int);
        posPos += myRaw(&rPos[posPos], computedRowMpz, sizeNum);
        
        ansPos.attr("class") = Rcpp::CharacterVector::create("bigz");
        return(ansPos);
    } else {
        if (computedRows > std::numeric_limits<int>::max())
            return Rcpp::wrap(computedRows);
        else
            return Rcpp::wrap(static_cast<int>(computedRows));
    }
}

void SetBounds(bool IsCount, SEXP Rlow, SEXP Rhigh,bool IsGmp, bool &bLower, 
               bool &bUpper, double &lower, double &upper, mpz_t *lowerMpz, mpz_t *upperMpz) {
    
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
}

void CheckBounds(bool IsGmp, double lower, double upper, double computedRows,
                 mpz_t lowerMpz, mpz_t upperMpz, mpz_t computedRowMpz) {
    
    if (IsGmp) {
        if (mpz_cmp(lowerMpz, computedRowMpz) >= 0 || mpz_cmp(upperMpz, computedRowMpz) > 0)
            Rcpp::stop("bounds cannot exceed the maximum number of possible results");
    } else {
        if (lower >= computedRows || upper > computedRows)
            Rcpp::stop("bounds cannot exceed the maximum number of possible results");
    }
}

void SetNumResults(bool IsGmp, bool bLower, bool bUpper, bool IsConstrained, bool &permNonTriv,
                   mpz_t *upperMpz, mpz_t *lowerMpz, double lower, double upper, double computedRows,
                   mpz_t &computedRowMpz, int &nRows, double &userNumRows) {
    
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
            permNonTriv = true;
            
            if (mpz_cmp_d(upperMpz[0], std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(upperMpz[0]);
        } else if (bLower) {
            mpz_sub(testBound, computedRowMpz, lowerMpz[0]);
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
            if (computedRows > std::numeric_limits<int>::max() && !IsConstrained)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = computedRows;
            
            if (!IsConstrained)
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

void SetRandomSampleMpz(SEXP RindexVec, SEXP RmySeed, std::size_t sampSize,
                        bool IsGmp, mpz_t &computedRowMpz, mpz_t *const myVec) {
    
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
                mpz_urandomm(myVec[i], seed_state, computedRowMpz);
            }
        }
        
        mpz_t maxGmp;
        mpz_init(maxGmp);
        mpz_set(maxGmp, myVec[0]);
        
        for (std::size_t i = 1; i < sampSize; ++i)
            if (mpz_cmp(myVec[i], maxGmp) > 0)
                mpz_set(maxGmp, myVec[i]);
            
        if (mpz_cmp(maxGmp, computedRowMpz) >= 0) {
            Rcpp::stop("One or more of the requested values in sampleVec "
                           "exceeds the maximum number of possible results");
        }
    }
}

// All functions below are exactly the same as the functions
// in CombPermUtils.cpp. The only difference is that they
// utilize the gmp library and deal mostly with mpz_t types

void NumPermsWithRepGmp(mpz_t result, const std::vector<int> &v) {
    
    mpz_set_ui(result, 1);
    std::vector<int> myLens = rleCpp(v);
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());
    
    const int myMax = myLens[0];
    const int numUni = myLens.size();
    
    for (int i = v.size(); i > myMax; --i)
        mpz_mul_ui(result, result, i);
    
    if (numUni > 1)
        for (int i = 1; i < numUni; ++i)
            for (int j = 2; j <= myLens[i]; ++j)
                mpz_divexact_ui(result, result, j);
}

void NumPermsNoRepGmp(mpz_t result, int n, int k) {
    
    mpz_set_ui(result, 1);
    
    for (int i = n, m = n - k; i > m; --i)
        mpz_mul_ui(result, result, i);
}

void nChooseKGmp(mpz_t result, int n, int k) {
    
    mpz_set_ui(result, 1);
    
    if (k != n && k != 0) {
        for (int i = (n - k + 1), d = 1; d <= k; ++i, ++d) {
            mpz_mul_ui(result, result, i);
            mpz_divexact_ui(result, result, d);
        }
    }
}

void NumCombsWithRepGmp(mpz_t result, int n, int r) {
    nChooseKGmp(result, n + r - 1, r);
}

bool onlyOneCombo(int n, int r, const std::vector<int> &Reps) {
    if (r < 1 || n <= 1)
        return true;
    
    if (r == n)
        if (std::accumulate(Reps.cbegin(), Reps.cend(), 0) == n)
            return true;
    
    return false;
}

void MultisetCombRowNumGmp(mpz_t result, int n, int r, 
                           const std::vector<int> &Reps) {
    
    if (!onlyOneCombo(n, r, Reps)) {
        const int r1 = r + 1;
        int myMax = r1;
        
        if (myMax > Reps[0] + 1)
            myMax = Reps[0] + 1;

        auto triangleVec = FromCpp14::make_unique<mpz_t[]>(r1);
        auto temp = FromCpp14::make_unique<mpz_t[]>(r1);
        
        for (int i = 0; i < r1; ++i) {
            mpz_init(triangleVec[i]);
            mpz_init(temp[i]);
        }
        
        for (int i = 0; i < myMax; ++i) {
            mpz_set_ui(triangleVec[i], 1);
            mpz_set_ui(temp[i], 1);
        }
        
        --myMax;
        int ind = 1;
        
        for (; myMax < r; ++ind) {
            int myMin = std::min(Reps[ind], r);
            
            for (int i = 1; i <= myMin; ++i)
                mpz_add(triangleVec[i], triangleVec[i], triangleVec[i - 1]);
            
            myMin = std::min(Reps[ind] + myMax, r);
            int j = 0;
            
            for (int i = (Reps[ind] + 1); i <= myMin; ++i, ++j) {
                mpz_add(triangleVec[i], triangleVec[i], triangleVec[i - 1]);
                mpz_sub(triangleVec[i], triangleVec[i], temp[j]);
                mpz_set(temp[j], triangleVec[j]);
            }
            
            for (; j <= myMin; ++j)
                mpz_set(temp[j], triangleVec[j]);
            
            myMax = myMin;
        }
        
        const int n1 = n - 1;
        mpz_t mySum, t;
        mpz_init(mySum);
        mpz_init(t);
        
        for (; ind < n1; ++ind) {
            mpz_set(t, triangleVec[r]);
            const int s = std::min(Reps[ind], r);
            
            for (int i = 1; i <= s; ++i)
                mpz_add(triangleVec[r], triangleVec[r], triangleVec[r - i]);
            
            mpz_set(mySum, triangleVec[r]);
            
            for (int i = r - 1; i >= s; --i) {
                mpz_sub(mySum, mySum, t);
                mpz_set(t, triangleVec[i]);
                mpz_add(mySum, mySum, triangleVec[i - s]);
                mpz_set(triangleVec[i], mySum);
            }
            
            for (int i = s - 1; i > 0; --i) {
                mpz_sub(mySum, mySum, t);
                mpz_set(t, triangleVec[i]);
                mpz_set(triangleVec[i], mySum);
            }
        }
        
        if (ind < n) {
            const int myMin2 = std::min(Reps[n1], r);
            
            for (int i = 1; i <= myMin2; ++i)
                mpz_add(triangleVec[r], triangleVec[r], triangleVec[r - i]);
        }
        
        mpz_set(result, triangleVec[r]);
        
        for (int i = 0; i < r1; ++i) {
            mpz_clear(triangleVec[i]);
            mpz_clear(temp[i]);
        }

        mpz_clear(mySum); mpz_clear(t);
    } else {
        mpz_set_ui(result, 1);
    }
}

void MultisetPermRowNumGmp(mpz_t result, int n, int r, const std::vector<int> &myReps) {
    
    const int sumFreqs = std::accumulate(myReps.cbegin(), myReps.cend(), 0);
    
    if (n < 2 || r < 1) {
        mpz_set_ui(result, 1);
    } else if (r > sumFreqs) {
        mpz_set_ui(result, 0);
    } else {
        const int n1 = n - 1;
        int maxFreq = *std::max_element(myReps.cbegin(), myReps.cend());
        
        std::vector<int> seqR(r);
        std::iota(seqR.begin(), seqR.end(), 1);
        
        mpz_t prodR;
        mpz_init(prodR);
        mpz_set_ui(prodR, 1);
        
        for (int i = 0; i < r; ++i)
            mpz_mul_ui(prodR, prodR, seqR[i]);
        
        const std::size_t uR1 = r + 1;
        const int myMax = (r < maxFreq) ? (r + 2) : (maxFreq + 2);
        auto cumProd = FromCpp14::make_unique<mpz_t[]>(myMax);
        auto resV = FromCpp14::make_unique<mpz_t[]>(uR1);
        
        for (int i = 0; i < myMax; ++i)
            mpz_init(cumProd[i]);
        
        // Equivalent to c(1, 1:myMax)
        mpz_set_ui(cumProd[0], 1);
        for (int i = 1; i < myMax; ++i)
            mpz_set_ui(cumProd[i], i);
        
        for (std::size_t i = 0; i < uR1; ++i) {
            mpz_init(resV[i]);
            mpz_set_ui(resV[i], 0);
        }
        
        for (int i = 1; i < myMax; ++i)
            mpz_mul(cumProd[i], cumProd[i], cumProd[i - 1]);
        
        int myMin = std::min(r, myReps[0]);
        
        for (int i = 0; i <= myMin; ++i)
            mpz_divexact(resV[i], prodR, cumProd[i]);
        
        mpz_t temp;
        mpz_init(temp);
        
        for (int i = 1; i < n1; ++i) {
            for (int j = r; j > 0; --j) {
                myMin = std::min(j, myReps[i]);
                mpz_set_ui(result, 0);
                
                for (int k = 0; k <= myMin; ++k) {
                    mpz_divexact(temp, resV[j - k], cumProd[k]);
                    mpz_add(result, result, temp);
                }
                
                mpz_set(resV[j], result);
            }
        }
        
        myMin = std::min(r, myReps[n1]);
        mpz_set_ui(result, 0);
        
        for (int k = 0; k <= myMin; ++k) {
            mpz_divexact(temp, resV[r - k], cumProd[k]);
            mpz_add(result, result, temp);
        }
        
        mpz_clear(temp); mpz_clear(prodR);
        
        for (int i = 0; i < myMax; ++i)
            mpz_clear(cumProd[i]);
        
        for (std::size_t i = 0; i < uR1; ++i)
            mpz_clear(resV[i]);
    }
}

void GetComputedRowMpz(mpz_t computedRowMpz, bool IsMultiset, bool IsComb, bool IsRep,
                       int n, int m, SEXP Rm, std::vector<int> &freqs, std::vector<int> &myReps) {
    
    if (IsMultiset) {
        if (IsComb) {
            MultisetCombRowNumGmp(computedRowMpz, n, m, myReps);
        } else {
            if (Rf_isNull(Rm) || m == static_cast<int>(freqs.size()))
                NumPermsWithRepGmp(computedRowMpz, freqs);
            else
                MultisetPermRowNumGmp(computedRowMpz, n, m, myReps);
        }
    } else {
        if (IsRep) {
            if (IsComb)
                NumCombsWithRepGmp(computedRowMpz, n, m);
            else
                mpz_ui_pow_ui(computedRowMpz, n, m);
        } else {
            if (IsComb)
                nChooseKGmp(computedRowMpz, n, m);
            else
                NumPermsNoRepGmp(computedRowMpz, n, m);
        }
    }
}
