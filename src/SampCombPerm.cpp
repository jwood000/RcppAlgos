#include "GmpDependUtils.h"
#include "RMatrix.h"
#include <RcppThread.h>

constexpr double dblDefault = 0;

template <typename typeRcpp, typename typeVector>
void SampleResults(const typeVector &v, std::size_t m, const std::vector<int> &myReps, 
                   std::size_t strtIdx, std::size_t endIdx, nthResultPtr nthResFun, 
                   const std::vector<double> &mySample, mpz_t *const myBigSamp, 
                   typeRcpp &sampleMatrix, int lenV, bool IsGmp) {
    
    if (IsGmp) {
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, dblDefault, myBigSamp[i], myReps);
            
            for (std::size_t j = 0; j < m; ++j)
                sampleMatrix(i, j) = v[z[j]];
        }
    } else {
        mpz_t mpzDefault;
        mpz_init(mpzDefault);
        
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[i], mpzDefault, myReps);
            
            for (std::size_t j = 0; j < m; ++j)
                sampleMatrix(i, j) = v[z[j]];
        }
        
        mpz_clear(mpzDefault);
    }
}

template <int RTYPE>
void SampNoThrdSafe(const Rcpp::Vector<RTYPE> &v, Rcpp::Matrix<RTYPE> &matRcpp,
                    const std::vector<int> &myReps, const std::vector<double> &mySample,
                    mpz_t *const myBigSamp, std::size_t m, std::size_t sampSize,
                    nthResultPtr nthResFun, int lenV, bool IsGmp, bool IsNamed) {
    
    SampleResults(v, m, myReps, 0, sampSize, nthResFun,
                  mySample, myBigSamp, matRcpp, lenV, IsGmp);
    
    if (IsNamed)
        SetSampleNames(IsGmp, sampSize, matRcpp, mySample, myBigSamp);
}

template <typename typeRcpp, typename typeVector>
void ParallelGlue(const std::vector<typeVector> &v, std::size_t m, const std::vector<int> &myReps, 
                  std::size_t strtIdx, std::size_t endIdx, nthResultPtr nthResFun, 
                  const std::vector<double> &mySample, mpz_t *const myBigSamp, 
                  typeRcpp &sampleMatrix, int lenV, bool IsGmp) {
    
    SampleResults(v, m, myReps, strtIdx, endIdx, nthResFun,
                  mySample, myBigSamp, sampleMatrix, lenV, IsGmp);
}

template <int RTYPE>
void SampleApplyFun(const Rcpp::Vector<RTYPE> &v, Rcpp::List &myList, std::size_t m,
                    const std::vector<int> &myReps, std::size_t sampSize,
                    const std::vector<double> &mySample, mpz_t *const myBigSamp,
                    SEXP func, SEXP rho, nthResultPtr nthResFun, bool IsNamed,
                    bool IsGmp, int lenV) {

    Rcpp::Vector<RTYPE> vectorPass(m);
    SEXP sexpFun = PROTECT(Rf_lang2(func, R_NilValue));
    
    if (IsGmp) {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, dblDefault, myBigSamp[i], myReps);
            
            for (std::size_t j = 0; j < m; ++j)
                vectorPass[j] = v[z[j]];
            
            SETCADR(sexpFun, vectorPass);
            myList[i] = Rf_eval(sexpFun, rho);
        }
    } else {
        mpz_t mpzDefault;
        mpz_init(mpzDefault);
        
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[i], mpzDefault, myReps);
            
            for (std::size_t j = 0; j < m; ++j)
                vectorPass[j] = v[z[j]];
            
            SETCADR(sexpFun, vectorPass);
            myList[i] = Rf_eval(sexpFun, rho);
        }
        
        mpz_clear(mpzDefault);
    }
    
    UNPROTECT(1);
    
    if (IsNamed)
        SetSampleNames(IsGmp, sampSize, myList, mySample, myBigSamp, false);
}

template <typename typeRcpp, typename typeElem>
void MasterSample(std::vector<typeElem> v, std::size_t m, const std::vector<int> &myReps,
                  std::size_t sampSize, nthResultPtr nthResFun, const std::vector<double> &mySample,
                  mpz_t *const myBigSamp, typeRcpp &matRcpp, int nThreads, bool Parallel,
                  bool IsNamed, bool IsGmp, int lenV) {
    
    if (Parallel) {
        RcppParallel::RMatrix<typeElem> parMat(matRcpp);
        RcppThread::ThreadPool pool(nThreads);
        int step = 0, stepSize = sampSize / nThreads;
        int nextStep = stepSize;

        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
            pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), 
                      std::cref(v), m, std::cref(myReps), step, nextStep, nthResFun, 
                      std::cref(mySample), myBigSamp, std::ref(parMat), lenV, IsGmp);
        }
        
        pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), 
                  std::cref(v), m, std::cref(myReps), step, sampSize, nthResFun, 
                  std::cref(mySample), myBigSamp, std::ref(parMat), lenV, IsGmp);
            
        pool.join();
    } else {
        SampleResults(v, m, myReps, 0, sampSize, nthResFun,
                      mySample, myBigSamp, matRcpp, lenV, IsGmp);
    }
    
    if (IsNamed)
        SetSampleNames(IsGmp, sampSize, matRcpp, mySample, myBigSamp);
}

// [[Rcpp::export]]
SEXP SampleRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, SEXP RFreqs, SEXP RindexVec, 
                bool IsComb, bool IsFactor, SEXP RmySeed, SEXP RNumSamp, 
                Rcpp::Function baseSample, SEXP stdFun, SEXP myEnv, SEXP Rparallel,
                SEXP RNumThreads, int maxThreads, SEXP RNamed) {
    
    int n, m = 0, lenFreqs = 0;
    bool IsMultiset = false;
    VecType myType = VecType::Integer;
    bool IsNamed = CleanConvert::convertLogical(RNamed, "namedSample");
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqs;
    
    bool IsRep = CleanConvert::convertLogical(Rrepetition, "repetition");
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    
    SetType(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    SetFreqsAndM(RFreqs, IsMultiset, myReps, IsRep, lenFreqs, freqs, Rm, n, m);
    
    const double computedRows = GetComputedRows(IsMultiset, IsComb, IsRep, n,
                                                m, Rm, lenFreqs, freqs, myReps);
    
    // sampleLimit defined in GmpDependUtils.h.. see comments in GmpCDependUtils.h
    bool IsGmp = computedRows > sampleLimit;
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);
    
    if (IsGmp) {
        GetComputedRowMpz(computedRowMpz, IsMultiset,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }
    
    std::size_t sampSize;
    std::vector<double> mySample;
    SetRandomSample(RindexVec, RNumSamp, sampSize,
                    IsGmp, computedRows, mySample, baseSample);
    
    const std::size_t bigSampSize = (IsGmp) ? sampSize : 1;
    auto myVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);
    
    for (std::size_t i = 0; i < bigSampSize; ++i)
        mpz_init(myVec[i]);
    
    SetRandomSampleMpz(RindexVec, RmySeed, sampSize, IsGmp, computedRowMpz, myVec.get());
    const bool applyFun = !Rf_isNull(stdFun) && !IsFactor;
    
    int nThreads = 1;
    const int limit = 2;
    SetThreads(Parallel, maxThreads, sampSize, myType, nThreads, RNumThreads, limit);
    
    Rcpp::XPtr<nthResultPtr> xpNth = putNthResPtrInXPtr(IsComb, IsMultiset, IsRep, IsGmp);
    const nthResultPtr nthResFun = *xpNth;
    
    if (applyFun) {
        if (!Rf_isFunction(stdFun))
            Rcpp::stop("FUN must be a function!");
        
        Rcpp::List myList(sampSize);
        
        switch (myType) {
            case VecType::Character : {
                Rcpp::CharacterVector charVec(Rcpp::clone(Rv));
                SampleApplyFun(charVec, myList, m, myReps, sampSize, mySample,
                               myVec.get(), stdFun, myEnv, nthResFun, IsNamed, IsGmp, n);
                break;
            } case VecType::Complex : {
                Rcpp::ComplexVector cmplxVec(Rcpp::clone(Rv));
                SampleApplyFun(cmplxVec, myList, m, myReps, sampSize, mySample,
                               myVec.get(), stdFun, myEnv, nthResFun, IsNamed, IsGmp, n);
                break;
            } case VecType::Raw : {
                Rcpp::RawVector rawVec(Rcpp::clone(Rv));
                SampleApplyFun(rawVec, myList, m, myReps, sampSize, mySample,
                               myVec.get(), stdFun, myEnv, nthResFun, IsNamed, IsGmp, n);
                break;
            } case VecType::Logical : {
                Rcpp::LogicalVector boolVec(Rcpp::clone(Rv));
                SampleApplyFun(boolVec, myList, m, myReps, sampSize, mySample,
                               myVec.get(), stdFun, myEnv, nthResFun, IsNamed, IsGmp, n);
                break;
            } case VecType::Integer : {
                Rcpp::IntegerVector intVec = Rcpp::wrap(vInt);
                SampleApplyFun(intVec, myList, m, myReps, sampSize, mySample,
                               myVec.get(), stdFun, myEnv, nthResFun, IsNamed, IsGmp, n);
                break;
            } default : {
                Rcpp::NumericVector numVec = Rcpp::wrap(vNum);
                SampleApplyFun(numVec, myList, m, myReps, sampSize, mySample,
                               myVec.get(), stdFun, myEnv, nthResFun, IsNamed, IsGmp, n);
                break;
            }
        }
        
        return myList;
    }
    
    switch (myType) {
        case VecType::Character : {
            Rcpp::CharacterVector charVec(Rcpp::clone(Rv));
            Rcpp::CharacterMatrix matChar = Rcpp::no_init_matrix(sampSize, m);
            SampNoThrdSafe(charVec, matChar, myReps, mySample, myVec.get(),
                           m, sampSize, nthResFun, n, IsGmp, IsNamed);
            return matChar;
        } case VecType::Complex : {
            Rcpp::ComplexVector cmplxVec(Rcpp::clone(Rv));
            Rcpp::ComplexMatrix matCmplx = Rcpp::no_init_matrix(sampSize, m);
            SampNoThrdSafe(cmplxVec, matCmplx, myReps, mySample, myVec.get(),
                           m, sampSize, nthResFun, n, IsGmp, IsNamed);
            return matCmplx;
        } case VecType::Raw : {
            Rcpp::RawVector rawVec(Rcpp::clone(Rv));
            Rcpp::RawMatrix matRaw = Rcpp::no_init_matrix(sampSize, m);
            SampNoThrdSafe(rawVec, matRaw, myReps, mySample, myVec.get(),
                           m, sampSize, nthResFun, n, IsGmp, IsNamed);
            return matRaw;
        } case VecType::Logical : {
            Rcpp::LogicalVector boolVec(Rcpp::clone(Rv));
            Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(sampSize, m);
            SampNoThrdSafe(boolVec, matBool, myReps, mySample, myVec.get(),
                           m, sampSize, nthResFun, n, IsGmp, IsNamed);
            return matBool;
        } case VecType::Integer : {
            Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(sampSize, m);
            MasterSample(vInt, m, myReps, sampSize, nthResFun, mySample, 
                         myVec.get(), matInt, nThreads, Parallel, IsNamed, IsGmp, n);
            if (IsFactor) {SetFactorClass(matInt, Rv);}
            return matInt;
        } default : {
            Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(sampSize, m);
            MasterSample(vNum, m, myReps, sampSize, nthResFun, mySample,
                         myVec.get(), matNum, nThreads, Parallel, IsNamed, IsGmp, n);
            return matNum;
        }
    }
}

