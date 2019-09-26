#include "NthResult.h"
#include "CleanConvert.h"
#include "CountGmp.h"
#include "RMatrix.h"
#include <RcppThread.h>

// Based off the internal limitations of sample, we
// cannot utilize the full range of 53-bit significand
// precision. Here is the condition from `do_dample2`:
//     if (!R_FINITE(dn) || dn < 0 || dn > 4.5e15 || (k > 0 && dn == 0))
// Here is the source (line ~1800): 
//     https://github.com/wch/r-source/blob/trunk/src/main/unique.c
const double sampleLimit = 4500000000000000.0;

template <typename typeRcpp, typename typeVector>
void SampleResults(const typeVector &v, std::size_t m, bool IsRep, const std::vector<int> &myReps,
                   std::size_t s, std::size_t n, bool IsGmp, bool IsComb,
                   const std::vector<double> &mySample, mpz_t *myBigSamp, typeRcpp &sampleMatrix) {

    const int lenV = v.size();
    std::vector<int> z(m);
    const bool IsMult = (static_cast<int>(myReps.size()) == lenV) ? true : false;
    
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
SEXP SampleApplyFun(const typeVector &v, std::size_t m, bool IsRep, bool IsGmp,
                    const std::vector<int> &myReps, std::size_t n, bool IsComb,
                    const std::vector<double> &mySample, mpz_t *myBigSamp, SEXP func, SEXP rho) {

    const int lenV = v.size();
    const bool IsMult = (static_cast<int>(myReps.size()) == lenV) ? true : false;
    std::vector<int> z(m);
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, n));
    SEXP sexpFun = PROTECT(Rf_lang2(func, R_NilValue));
    typeVector vectorPass(m);
    
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
                Rcpp::Function baseSample, SEXP stdFun, SEXP myEnv, 
                SEXP Rparallel, SEXP RNumThreads, int maxThreads) {
    
    int n, m = 0, lenFreqs = 0;
    bool IsMultiset, IsInteger, IsCharacter, IsLogical;
    IsCharacter = IsInteger = IsLogical = false;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqsExpanded;
    Rcpp::CharacterVector rcppChar;
    
    bool IsRepetition = CleanConvert::convertLogical(Rrepetition, "repetition");
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    
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
    }
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        myReps.push_back(1);
    } else {
        IsRepetition = false;
        CleanConvert::convertVector(RFreqs, myReps, "freqs");
        int testTrivial = std::accumulate(myReps.cbegin(), myReps.cend(), 0);
        lenFreqs = static_cast<int>(myReps.size());
        
        if (testTrivial > lenFreqs) {
            IsMultiset = true;
            
            for (int i = 0; i < lenFreqs; ++i)
                for (int j = 0; j < myReps[i]; ++j)
                    freqsExpanded.push_back(i);
        } else {
            IsMultiset = false;
            freqsExpanded = myReps;
        }
    }
    
    if (Rf_isNull(Rm)) {
        if (!freqsExpanded.empty()) {
            m = freqsExpanded.size();
        } else {
            Rcpp::stop("m and freqs cannot both be NULL");
        }
    } else {
        if (Rf_length(Rm) > 1)
            Rcpp::stop("length of m must be 1");
        
        CleanConvert::convertPrimitive(Rm, m, "m");
    }
    
    SetValues(IsCharacter, IsLogical, IsInteger, rcppChar, vInt, vNum, n, Rv);
    
    if (IsFactor)
        IsLogical = IsCharacter = IsInteger = false;
    
    const double computedRows = GetComputedRows(IsMultiset, IsComb, IsRepetition, n,
                                                m, Rm, lenFreqs, freqsExpanded, myReps);
    
    // sampleLimit defined as const above... see comments for more details
    bool IsGmp = computedRows > sampleLimit;
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);
    
    if (IsGmp) {
        GetComputedRowMpz(computedRowMpz, IsMultiset, IsComb, 
                          IsRepetition, n, m, Rm, freqsExpanded, myReps);
    }
    
    std::size_t sampSize;
    std::vector<double> mySample;
    SetRandomSample(RindexVec, RNumSamp, sampSize, IsGmp, computedRows, mySample, baseSample);
    
    std::size_t gmpSize = (IsGmp) ? sampSize : 1;
    auto myVec = FromCpp14::make_unique<mpz_t[]>(gmpSize);
    SetRandomSampleMpz(RindexVec, RmySeed, sampSize, IsGmp, computedRowMpz, myVec.get());
    
    bool applyFun = !Rf_isNull(stdFun) && !IsFactor;
    int nThreads = 1;
    const int limit = 2;
    SetThread(Parallel, maxThreads, sampSize, IsCharacter, nThreads, RNumThreads, limit);
    
    if (applyFun) {
        if (!Rf_isFunction(stdFun))
            Rcpp::stop("FUN must be a function!");
        
        if (IsCharacter) {
            return SampleApplyFun(rcppChar, m, IsRepetition, IsGmp, myReps, 
                                  sampSize, IsComb, mySample, myVec.get(), stdFun, myEnv);
        } else if (IsLogical || IsInteger) {
            Rcpp::IntegerVector rcppVInt(vInt.begin(), vInt.end());
            return SampleApplyFun(rcppVInt, m, IsRepetition, IsGmp, myReps, 
                                  sampSize, IsComb, mySample, myVec.get(), stdFun, myEnv);
        } else {
            Rcpp::NumericVector rcppVNum(vNum.begin(), vNum.end());
            return SampleApplyFun(rcppVNum, m, IsRepetition, IsGmp, myReps, 
                                  sampSize, IsComb, mySample, myVec.get(), stdFun, myEnv);
        }
    }
    
    if (Parallel) {
        RcppThread::ThreadPool pool(nThreads);
        int step = 0, stepSize = sampSize / nThreads;
        int nextStep = stepSize;
        
        if (IsLogical) {
            Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(sampSize, m);
            RcppParallel::RMatrix<int> parBool(matBool);

            for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                pool.push(std::cref(SampleResults<RcppParallel::RMatrix<int>, std::vector<int>>), vInt, m,
                          IsRepetition, myReps, step, nextStep, IsGmp, IsComb, mySample, myVec.get(), std::ref(parBool));
            }
            
            pool.push(std::cref(SampleResults<RcppParallel::RMatrix<int>, std::vector<int>>), vInt, m,
                      IsRepetition, myReps, step, sampSize, IsGmp, IsComb, mySample, myVec.get(), std::ref(parBool));
            
            pool.join();
            return matBool;
            
        } else if (IsFactor || IsInteger) {
            Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(sampSize, m);
            RcppParallel::RMatrix<int> parInt(matInt);
            
            for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                pool.push(std::cref(SampleResults<RcppParallel::RMatrix<int>, std::vector<int>>), vInt, m,
                          IsRepetition, myReps, step, nextStep, IsGmp, IsComb, mySample, myVec.get(), std::ref(parInt));
            }
            
            pool.push(std::cref(SampleResults<RcppParallel::RMatrix<int>, std::vector<int>>), vInt, m,
                      IsRepetition, myReps, step, sampSize, IsGmp, IsComb, mySample, myVec.get(), std::ref(parInt));
            
            pool.join();
            
            if (IsFactor) {
                Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
                Rcpp::CharacterVector myClass = testFactor.attr("class");
                Rcpp::CharacterVector myLevels = testFactor.attr("levels");
                matInt.attr("class") = myClass;
                matInt.attr("levels") = myLevels;
            }
            
            return matInt;
            
        } else {
            Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(sampSize, m);
            RcppParallel::RMatrix<double> parNum(matNum);
            
            for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                pool.push(std::cref(SampleResults<RcppParallel::RMatrix<double>, std::vector<double>>), vNum, m,
                          IsRepetition, myReps, step, nextStep, IsGmp, IsComb, mySample, myVec.get(), std::ref(parNum));
            }
            
            pool.push(std::cref(SampleResults<RcppParallel::RMatrix<double>, std::vector<double>>), vNum, m,
                      IsRepetition, myReps, step, sampSize, IsGmp, IsComb, mySample, myVec.get(), std::ref(parNum));
            
            pool.join();
            return matNum;
        }
    }
    
    if (IsCharacter) {
        Rcpp::CharacterMatrix matChar = Rcpp::no_init_matrix(sampSize, m);
        SampleResults(rcppChar, m, IsRepetition, myReps, 0, sampSize, IsGmp, IsComb, mySample, myVec.get(), matChar);
        return matChar;
    } else if (IsLogical) {
        Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(sampSize, m);
        SampleResults(vInt, m, IsRepetition, myReps, 0, sampSize, IsGmp, IsComb, mySample, myVec.get(), matBool);
        return matBool;
    } else if (IsFactor || IsInteger) {
        Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(sampSize, m);
        SampleResults(vInt, m, IsRepetition, myReps, 0, sampSize, IsGmp, IsComb, mySample, myVec.get(), matInt);

        if (IsFactor) {
            Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
            Rcpp::CharacterVector myClass = testFactor.attr("class");
            Rcpp::CharacterVector myLevels = testFactor.attr("levels");
            matInt.attr("class") = myClass;
            matInt.attr("levels") = myLevels;
        }
        
        return matInt;
    } else {
        Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(sampSize, m);
        SampleResults(vNum, m, IsRepetition, myReps, 0, sampSize, IsGmp, IsComb, mySample, myVec.get(), matNum);
        return matNum;
    }
}

