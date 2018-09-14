#include <NthResult.h>
#include <CleanConvert.h>
#include <importExportMPZ.h>
#include <CombPermUtils.h>
#include <CountGmp.h>

// Based off the internal limitations of sample, we
// cannot utilize the full range of 53-bit significand
// precision. Here is the condition from `do_dample2`:
//     if (!R_FINITE(dn) || dn < 0 || dn > 4.5e15 || (k > 0 && dn == 0))
// Here is the source (line 1791): 
//     https://github.com/wch/r-source/blob/trunk/src/main/unique.c

static gmp_randstate_t seed_state;
static int seed_init = 0;
const double sampleLimit = 4500000000000000.0;

template <typename typeRcpp, typename typeVector>
void SampleResults(typeVector v, unsigned long int m, bool IsRep, std::vector<int> myReps,
                   unsigned long int s, unsigned long int n, bool IsGmp, bool IsComb,
                   std::vector<double> mySample, mpz_t myBigSamp[], typeRcpp sampleMatrix) {

    int lenV = v.size();
    std::vector<int> z(m);
    bool IsMult = false;
    
    if ((int) myReps.size() == lenV)
        IsMult = true;
    
    if (IsGmp) {
        if (IsComb) {
            for (std::size_t i = s; i < n; ++i) {
                z = nthCombinationGmp(lenV, m, myBigSamp[i], IsRep, IsMult, myReps);
                for (std::size_t j = 0; j < m; ++j)
                    sampleMatrix(i, j) = v[z[j]];
            }
        } else {
            for (std::size_t i = s; i < n; ++i) {
                z = nthPermutationGmp(lenV, m, myBigSamp[i], IsRep, IsMult, myReps, myReps);
                for (std::size_t j = 0; j < m; ++j)
                    sampleMatrix(i, j) = v[z[j]];
            }
        }
    } else {
        if (IsComb) {
            for (std::size_t i = s; i < n; ++i) {
                z = nthCombination(lenV, m, mySample[i] - 1, IsRep, IsMult, myReps);
                for (std::size_t j = 0; j < m; ++j)
                    sampleMatrix(i, j) = v[z[j]];
            }
        } else {
            for (std::size_t i = s; i < n; ++i) {
                z = nthPermutation(lenV, m, mySample[i] - 1, IsRep, IsMult, myReps, myReps);
                for (std::size_t j = 0; j < m; ++j)
                    sampleMatrix(i, j) = v[z[j]];
            }
        }
    }
}

template <typename typeVector>
SEXP SampleApplyFun(typeVector &v, unsigned long int m, bool IsRep, bool IsGmp,
                    std::vector<int> myReps, unsigned long int n, bool IsComb,
                    std::vector<double> &mySample, mpz_t myBigSamp[], SEXP func, SEXP rho) {

    int lenV = v.size();
    std::vector<int> z(m);
    bool IsMult = false;
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, n));
    SEXP sexpFun = PROTECT(Rf_lang2(func, R_NilValue));
    typeVector vectorPass(m);
    
    if ((int) myReps.size() == lenV)
        IsMult = true;
    
    if (IsGmp) {
        if (IsComb) {
            for (std::size_t i = 0; i < n; ++i) {
                z = nthCombinationGmp(lenV, m, myBigSamp[i], IsRep, IsMult, myReps);
                for (std::size_t j = 0; j < m; ++j)
                    vectorPass[j] = v[z[j]];
                
                SETCADR(sexpFun, vectorPass);
                SET_VECTOR_ELT(ans, i, Rf_eval(sexpFun, rho));
            }
        } else {
            for (std::size_t i = 0; i < n; ++i) {
                z = nthPermutationGmp(lenV, m, myBigSamp[i], IsRep, IsMult, myReps, myReps);
                for (std::size_t j = 0; j < m; ++j)
                    vectorPass[j] = v[z[j]];
                
                SETCADR(sexpFun, vectorPass);
                SET_VECTOR_ELT(ans, i, Rf_eval(sexpFun, rho));
            }
        }
    } else {
        if (IsComb) {
            for (std::size_t i = 0; i < n; ++i) {
                z = nthCombination(lenV, m, mySample[i] - 1, IsRep, IsMult, myReps);
                for (std::size_t j = 0; j < m; ++j)
                    vectorPass[j] = v[z[j]];
                
                SETCADR(sexpFun, vectorPass);
                SET_VECTOR_ELT(ans, i, Rf_eval(sexpFun, rho));
            }
        } else {
            for (std::size_t i = 0; i < n; ++i) {
                z = nthPermutation(lenV, m, mySample[i] - 1, IsRep, IsMult, myReps, myReps);
                for (std::size_t j = 0; j < m; ++j)
                    vectorPass[j] = v[z[j]];
                
                SETCADR(sexpFun, vectorPass);
                SET_VECTOR_ELT(ans, i, Rf_eval(sexpFun, rho));
            }
        }
    }
    
    UNPROTECT(2);
    return ans;
}

// [[Rcpp::export]]
SEXP SampleRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, SEXP RFreqs, SEXP RindexVec,
                bool IsComb, bool IsFactor, SEXP RmySeed, SEXP RNumSamp, 
                Rcpp::Function baseSample, SEXP stdFun, SEXP myEnv) {
    
    int n, m1, m2, m = 0, lenFreqs = 0;
    bool IsRepetition, IsLogical;
    bool IsMultiset, IsInteger, IsCharacter;
    IsCharacter = IsInteger = IsLogical = false;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqsExpanded;
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
            IsLogical = IsInteger = false;
            break;
        }
    }
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        myReps.push_back(1);
    } else {
        IsMultiset = true;
        CleanConvert::convertVector(RFreqs, myReps, "freqs must be of type numeric or integer");
        
        lenFreqs = (int) myReps.size();
        for (int i = 0; i < lenFreqs; ++i) {
            if (myReps[i] < 1) 
                Rcpp::stop("each element in freqs must be a positive number");
            
            for (int j = 0; j < myReps[i]; ++j)
                freqsExpanded.push_back(i);
        }
    }
    
    if (Rf_isNull(Rm)) {
        if (IsMultiset) {
            m = freqsExpanded.size();
        } else {
            Rcpp::stop("m and freqs cannot both be NULL");
        }
    } else {
        if (Rf_length(Rm) > 1)
            Rcpp::stop("length of m must be 1");
        
        CleanConvert::convertPrimitive(Rm, m, "m must be of type numeric or integer");
    }
    
    if (m < 1)
        Rcpp::stop("m must be positive");
    
    IsRepetition = Rcpp::as<bool>(Rrepetition);
    double seqEnd;
    
    if (IsCharacter) {
        rcppChar = Rcpp::as<Rcpp::CharacterVector>(Rv);
        n = rcppChar.size();
    } else if (IsLogical) {
        vInt = Rcpp::as<std::vector<int> >(Rv);
        n = vInt.size();
    } else {
        if (Rf_length(Rv) == 1) {
            seqEnd = Rcpp::as<double>(Rv);
            if (Rcpp::NumericVector::is_na(seqEnd)) {seqEnd = 1;}
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            Rcpp::IntegerVector vTemp = Rcpp::seq(m1, m2);
            IsInteger = true;
            vNum = Rcpp::as<std::vector<double> >(vTemp);
        } else {
            vNum = Rcpp::as<std::vector<double> >(Rv);
        }
        
        n = vNum.size();
    }
    
    if (IsInteger) {
        for (int i = 0; i < n && IsInteger; ++i)
            if (Rcpp::NumericVector::is_na(vNum[i]))
                IsInteger = false;
            
        if (IsInteger)
            vInt.assign(vNum.begin(), vNum.end());
    }
    
    if (IsFactor)
        IsLogical = IsCharacter = IsInteger = false;
    
    double computedRows;
    
    if (IsMultiset) {
        if (n != lenFreqs)
            Rcpp::stop("the length of freqs must equal the length of v");
        
        if (m > (int) freqsExpanded.size())
            m = freqsExpanded.size();
        
        if (IsComb) {
            computedRows = MultisetCombRowNum(n, m, myReps);
        } else {
            if (Rf_isNull(Rm) || m == (int) freqsExpanded.size())
                computedRows = NumPermsWithRep(freqsExpanded);
            else
                computedRows = MultisetPermRowNum(n, m, myReps);
        }
    } else {
        if (IsRepetition) {
            if (IsComb)
                computedRows = NumCombsWithRep(n, m);
            else
                computedRows = std::pow((double) n, (double) m);
        } else {
            if (m > n)
                Rcpp::stop("m must be less than or equal to the length of v");
            
            if (IsComb)
                computedRows = nChooseK(n, m);
            else
                computedRows = NumPermsNoRep(n, m);
        }
    }
    
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);
    bool IsGmp = false;
    
    // sampleLimit defined as const above... see comments for more details
    if (computedRows > sampleLimit) {
        IsGmp = true;
        if (IsMultiset) {
            if (IsComb) {
                MultisetCombRowNumGmp(computedRowMpz, n, m, myReps);
            } else {
                if (Rf_isNull(Rm) || m == (int) freqsExpanded.size())
                    NumPermsWithRepGmp(computedRowMpz, freqsExpanded);
                else
                    MultisetPermRowNumGmp(computedRowMpz, n, m, myReps);
            }
        } else {
            if (IsRepetition) {
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
    
    unsigned long int sampSize;
    std::vector<double> mySample;
    
    if (Rf_isNull(RindexVec)) {
        if (Rf_isNull(RNumSamp))
            Rcpp::stop("n and sampleVec cannot both be NULL");
    
        if (!IsGmp) {
            if (!Rf_isNumber(RNumSamp))
                Rcpp::stop("n must be a number");
            
            if (Rf_length(RNumSamp) > 1)
                Rcpp::stop("length of n must be 1. For specific combinations, use sampleVec.");
            
            double nPass;
            CleanConvert::convertPrimitive(RNumSamp, nPass, "n must be a number");
            
            if (nPass > computedRows)
                Rcpp::stop("n exceeds the maximum number of possible results");
                
            if (nPass > INT_MAX)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1");
            
            Rcpp::NumericVector tempSamp = baseSample(computedRows, nPass);
            mySample = Rcpp::as<std::vector<double> >(tempSamp);
        }
    } else if (!IsGmp) {
        CleanConvert::convertVector(RindexVec, mySample, 
                                    "sampleVec must be convertible to a real number");
        if (mySample.size() == 1)
            if (mySample[0] == 0)
                Rcpp::stop("sampleVec must be convertible to a real number");
    }
    
    // Initialize myVec to silence compiler
    mpz_t *myVec;
    myVec = (mpz_t *) malloc(sizeof(mpz_t));
    mpz_init(myVec[0]);
    
    if (IsGmp) {
        mpz_t maxGmp;
        mpz_init(maxGmp);
        
        if (!Rf_isNull(RindexVec)) {
            switch (TYPEOF(RindexVec)) {
                case RAWSXP: {
                    const char* raw = (char*)RAW(RindexVec);
                    sampSize = ((int*)raw)[0];
                    break;
                }
                default:
                    sampSize = LENGTH(RindexVec);
            }
            
            mpz_clear(myVec[0]);
            myVec = (mpz_t *) malloc(sizeof(mpz_t) * sampSize);
            for (std::size_t i = 0; i < sampSize; ++i)
                mpz_init(myVec[i]);
            
            createMPZArray(RindexVec, myVec, sampSize);

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
                createMPZArray(RmySeed, mpzSeed, 1);
                gmp_randseed(seed_state, mpzSeed[0]);
                mpz_clear(mpzSeed[0]);
            }
            
            double dblVSize;
            CleanConvert::convertPrimitive(RNumSamp, dblVSize, "n must be convertible to a number");
            
            if (dblVSize > INT_MAX)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1");
            
            if (dblVSize < 1)
                Rcpp::stop("n must be a positive number");
            
            sampSize = (unsigned long int) dblVSize;
            myVec = (mpz_t *) malloc(sizeof(mpz_t) * sampSize);
            
            // random number is between 0 and gmpRows[0] - 1
            // so we need to add 1 to each element
            for (std::size_t i = 0; i < sampSize; ++i) {
                mpz_init(myVec[i]);
                mpz_urandomm(myVec[i], seed_state, computedRowMpz);
            }
        }
        
        mpz_set(maxGmp, myVec[0]);
        
        for (std::size_t i = 0; i < sampSize; ++i)
            if (mpz_cmp(myVec[i], maxGmp) > 0)
                mpz_set(maxGmp, myVec[i]);
        
        if (mpz_cmp(maxGmp, computedRowMpz) >= 0)
            Rcpp::stop("One or more of the requested values in sampleVec "
                           "exceeds the maximum number of possible results");
    } else {
        sampSize = mySample.size();
        if (sampSize > INT_MAX)
            Rcpp::stop("The number of rows cannot exceed 2^31 - 1");
        
        double myMax = *std::max_element(mySample.begin(), mySample.end());
        if (myMax > computedRows)
            Rcpp::stop("One or more of the requested values in sampleVec "
                           "exceeds the maximum number of possible results");
    }
    
    bool applyFun = !Rf_isNull(stdFun) && !IsFactor;
    
    if (applyFun) {
        if (!Rf_isFunction(stdFun))
            Rcpp::stop("FUN must be a function!");
        
        if (IsCharacter) {
            return SampleApplyFun(rcppChar, m, IsRepetition, IsGmp, myReps, 
                                  sampSize, IsComb, mySample, myVec, stdFun, myEnv);
        } else if (IsLogical || IsInteger) {
            Rcpp::IntegerVector rcppVInt(vInt.begin(), vInt.end());
            return SampleApplyFun(rcppVInt, m, IsRepetition, IsGmp, myReps, 
                                  sampSize, IsComb, mySample, myVec, stdFun, myEnv);
        } else {
            Rcpp::NumericVector rcppVNum(vNum.begin(), vNum.end());
            return SampleApplyFun(rcppVNum, m, IsRepetition, IsGmp, myReps, 
                                  sampSize, IsComb, mySample, myVec, stdFun, myEnv);
        }
    }
    
    if (IsCharacter) {
        Rcpp::CharacterMatrix matChar = Rcpp::no_init_matrix(sampSize, m);
        SampleResults(rcppChar, m, IsRepetition, myReps, 0, sampSize,
                      IsGmp, IsComb, mySample, myVec, matChar);
        return matChar;
    } else if (IsLogical) {
        Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(sampSize, m);
        SampleResults(vInt, m, IsRepetition, myReps, 0, sampSize,
                      IsGmp, IsComb, mySample, myVec, matBool);
        return matBool;
    } else if (IsFactor) {
        Rcpp::IntegerMatrix factorMat = Rcpp::no_init_matrix(sampSize, m);
        Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
        Rcpp::CharacterVector myClass = testFactor.attr("class");
        Rcpp::CharacterVector myLevels = testFactor.attr("levels");

        SampleResults(vInt, m, IsRepetition, myReps, 0, sampSize,
                      IsGmp, IsComb, mySample, myVec, factorMat);

        factorMat.attr("class") = myClass;
        factorMat.attr("levels") = myLevels;

        return factorMat;
    } else if (IsInteger) {
        Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(sampSize, m);
        SampleResults(vInt, m, IsRepetition, myReps, 0, sampSize,
                      IsGmp, IsComb, mySample, myVec, matInt);
        return matInt;
    } else {
        Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(sampSize, m);
        SampleResults(vNum, m, IsRepetition, myReps, 0, sampSize,
                      IsGmp, IsComb, mySample, myVec, matNum);
        return matNum;
    }
}

