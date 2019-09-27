#include "NthResult.h"
#include "CleanConvert.h"
#include "CountGmp.h"
#include "RMatrix.h"
#include <RcppThread.h>

template <typename typeRcpp, typename typeVector>
void ParallelGlue(std::vector<typeVector> v, std::size_t m, const std::vector<int> &myReps, 
                  std::size_t strtIdx, std::size_t endIdx, nthResutlPtr nthResFun, 
                  const std::vector<double> &mySample, mpz_t *const myBigSamp, typeRcpp &sampleMatrix) {
    SampleResults(v, m, myReps, strtIdx, endIdx, nthResFun, mySample, myBigSamp, sampleMatrix);
}

template <typename typeRcpp, typename typeVector>
void SampleResults(const typeVector &v, std::size_t m, const std::vector<int> &myReps, 
                   std::size_t strtIdx, std::size_t endIdx, nthResutlPtr nthResFun, 
                   const std::vector<double> &mySample, mpz_t *const myBigSamp, typeRcpp &sampleMatrix) {

    const int lenV = v.size();
    std::vector<int> z(m);
    
    for (std::size_t i = strtIdx; i < endIdx; ++i) {
        z = nthResFun(lenV, m, mySample[i], myBigSamp[i], myReps);
        for (std::size_t j = 0; j < m; ++j)
            sampleMatrix(i, j) = v[z[j]];
    }
}

template <typename typeVector>
SEXP SampleApplyFun(const typeVector &v, std::size_t m, const std::vector<int> &myReps,
                    std::size_t sampSize, const std::vector<double> &mySample, 
                    mpz_t *const myBigSamp, SEXP func, SEXP rho, nthResutlPtr nthResFun) {

    const int lenV = v.size();
    std::vector<int> z(m);
    typeVector vectorPass(m);
    
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, sampSize));
    SEXP sexpFun = PROTECT(Rf_lang2(func, R_NilValue));
    
    for (std::size_t i = 0; i < sampSize; ++i) {
        z = nthResFun(lenV, m, mySample[i], myBigSamp[i], myReps);
        for (std::size_t j = 0; j < m; ++j)
            vectorPass[j] = v[z[j]];
        
        SETCADR(sexpFun, vectorPass);
        SET_VECTOR_ELT(ans, i, Rf_eval(sexpFun, rho));
    }

    UNPROTECT(2);
    return ans;
}

template <typename typeRcpp, typename typeElem>
void MasterSample(std::vector<typeElem> v, std::size_t m, const std::vector<int> &myReps,
                  std::size_t sampSize, nthResutlPtr nthResFun, const std::vector<double> &mySample,
                  mpz_t *const myBigSamp, typeRcpp &matRcpp, int nThreads, bool Parallel) {
    
    if (Parallel) {
        RcppParallel::RMatrix<typeElem> parMat(matRcpp);
        RcppThread::ThreadPool pool(nThreads);
        int step = 0, stepSize = sampSize / nThreads;
        int nextStep = stepSize;

        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
            pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), v,
                      m, myReps, step, nextStep, nthResFun, mySample, myBigSamp, std::ref(parMat));
        }
        
        pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), v,
                  m, myReps, step, sampSize, nthResFun, mySample, myBigSamp, std::ref(parMat));
            
        pool.join();
    } else {
        SampleResults(v, m, myReps, 0, sampSize, nthResFun, mySample, myBigSamp, matRcpp);
    }
}

// [[Rcpp::export]]
SEXP SampleRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, SEXP RFreqs, SEXP RindexVec, bool IsComb,
                bool IsFactor, SEXP RmySeed, SEXP RNumSamp, Rcpp::Function baseSample, SEXP stdFun,
                SEXP myEnv, SEXP Rparallel, SEXP RNumThreads, int maxThreads) {
    
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
    
    // sampleLimit defined in CountGmp.h.. see comments in CountGmp.h for more details
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
    auto myVec = FromCpp14::make_unique<mpz_t[]>(sampSize);
    
    for (std::size_t i = 0; i < sampSize; ++i)
        mpz_init(myVec[i]);
    
    SetRandomSampleMpz(RindexVec, RmySeed, sampSize, IsGmp, computedRowMpz, myVec.get());
    bool applyFun = !Rf_isNull(stdFun) && !IsFactor;
    int nThreads = 1;
    const int limit = 2;
    
    SetThreads(Parallel, maxThreads, sampSize, IsCharacter, nThreads, RNumThreads, limit);
    Rcpp::XPtr<nthResutlPtr> xpNth = putNthResPtrInXPtr(IsMultiset, IsRepetition, IsGmp, IsComb);
    nthResutlPtr nthResFun = *xpNth;
    
    if (applyFun) {
        if (!Rf_isFunction(stdFun))
            Rcpp::stop("FUN must be a function!");

        if (IsCharacter) {
            return SampleApplyFun(rcppChar, m, myReps, sampSize, 
                                  mySample, myVec.get(), stdFun, myEnv, nthResFun);
        } else if (IsLogical || IsInteger) {
            Rcpp::IntegerVector rcppVInt(vInt.begin(), vInt.end());
            return SampleApplyFun(rcppVInt, m, myReps, sampSize, 
                                  mySample, myVec.get(), stdFun, myEnv, nthResFun);
        } else {
            Rcpp::NumericVector rcppVNum(vNum.begin(), vNum.end());
            return SampleApplyFun(rcppVNum, m, myReps, sampSize, 
                                  mySample, myVec.get(), stdFun, myEnv, nthResFun);
        }
    }
    
    if (IsCharacter) {
        Rcpp::CharacterMatrix matChar = Rcpp::no_init_matrix(sampSize, m);
        SampleResults(rcppChar, m, myReps, 0, sampSize, nthResFun, mySample, myVec.get(), matChar);
        return matChar;
    } else if (IsLogical) {
        Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(sampSize, m);
        MasterSample(vInt, m, myReps, sampSize, nthResFun, 
                     mySample, myVec.get(), matBool, nThreads, Parallel);
        return matBool;
    } else if (IsFactor || IsInteger) {
        Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(sampSize, m);
        MasterSample(vInt, m, myReps, sampSize, nthResFun, 
                     mySample, myVec.get(), matInt, nThreads, Parallel);

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
        MasterSample(vNum, m, myReps, sampSize, nthResFun, 
                     mySample, myVec.get(), matNum, nThreads, Parallel);
        return matNum;
    }
}

