#include "Cpp14MakeUnique.h"
#include "SampleCombPerm.h"
#include "SampleApplyFun.h"
#include "ComputedCount.h"
#include "SampleResults.h"
#include <thread>

template <typename T>
void SampNoThrdSafe(T* sampleMatrix, SEXP res, const std::vector<T> &v,
                    const std::vector<double> &mySample,
                    mpz_t *const myBigSamp, const std::vector<int> &myReps,
                    nthResultPtr nthResFun, int m, int sampSize,
                    int lenV, bool IsGmp, bool IsNamed) {

    SampleResults(sampleMatrix, v, mySample, myBigSamp, myReps,
                  nthResFun, m, sampSize, lenV, IsGmp);

    if (IsNamed) {
        SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp);
    }
}

template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T> &sampleMatrix,
                  const std::vector<T> &v,
                  const std::vector<double> &mySample,
                  mpz_t *const myBigSamp, const std::vector<int> &myReps,
                  nthResultPtr nthResFun, int m, int strtIdx, int endIdx,
                  int lenV, bool IsGmp) {

    SampleResults(sampleMatrix, v, mySample, myBigSamp, myReps,
                  nthResFun, m, strtIdx, endIdx, lenV, IsGmp);
}

template <typename T>
void ThreadSafeSample(T* mat, SEXP res, const std::vector<T> &v,
                      const std::vector<double> &mySample,
                      mpz_t *const myBigSamp, const std::vector<int> &myReps,
                      nthResultPtr nthResFun, int m, int sampSize,
                      int nThreads, bool Parallel, bool IsNamed,
                      bool IsGmp, int lenV) {

    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(mat, sampSize, m);
        std::vector<std::thread> threads;

        int step = 0;
        int stepSize = sampSize / nThreads;
        int nextStep = stepSize;

        for (int j = 0; j < (nThreads - 1);
             ++j, step += stepSize, nextStep += stepSize) {

            threads.emplace_back(std::cref(ParallelGlue<T>),
                                 std::ref(parMat), std::cref(v),
                                 std::cref(mySample), myBigSamp,
                                 std::cref(myReps), nthResFun,
                                 m, step, nextStep, lenV, IsGmp);
        }

        threads.emplace_back(std::cref(ParallelGlue<T>), std::ref(parMat),
                             std::cref(v), std::cref(mySample), myBigSamp,
                             std::cref(myReps), nthResFun, m, step, sampSize,
                             lenV, IsGmp);

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        SampleResults(mat, v, mySample, myBigSamp, myReps,
                      nthResFun, m, sampSize, lenV, IsGmp);
    }

    if (IsNamed) {
        SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp);
    }
}

SEXP SampleCombPerm(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                    SEXP RindexVec, SEXP RIsComb, SEXP RmySeed,
                    SEXP RNumSamp, SEXP baseSample, SEXP stdFun, SEXP myEnv,
                    SEXP Rparallel, SEXP RNumThreads, SEXP RmaxThreads,
                    SEXP RNamed, SEXP RFunVal) {

    int n = 0;
    int m = 0;
    int nThreads = 1;
    int maxThreads = 1;

    VecType myType = VecType::Integer;
    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");
    bool IsNamed = CleanConvert::convertLogical(RNamed, "namedSample");

    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;
    std::vector<double> vNum;

    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    const bool IsComb = CleanConvert::convertLogical(RIsComb, "IsComb");
    bool IsMult = false;

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum,
              Rv, RFreqs, Rm, n, m, IsMult, IsRep);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp = (computedRows > SampleLimit);

    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult, IsComb,
                          IsRep, n, m, Rm, freqs, myReps);
    }

    int sampSize;
    std::vector<double> mySample;
    SetRandomSample(RindexVec, RNumSamp, sampSize, IsGmp,
                    computedRows, mySample, baseSample, myEnv);

    const int bigSampSize = (IsGmp) ? sampSize : 1;
    auto myVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);

    for (int i = 0; i < bigSampSize; ++i) {
        mpz_init(myVec[i]);
    }

    SetRandomSampleMpz(RindexVec, RmySeed, sampSize,
                       IsGmp, computedRowsMpz, myVec.get());

    const int limit = 2;
    SetThreads(Parallel, maxThreads, sampSize,
               myType, nThreads, RNumThreads, limit);

    const nthResultPtr nthResFun = GetNthResultFunc(IsComb, IsMult,
                                                    IsRep, IsGmp);

    if (!Rf_isNull(stdFun) && !Rf_isFactor(Rv)) {
        if (!Rf_isFunction(stdFun)) {
            Rf_error("FUN must be a function!");
        }

        switch (myType) {
            case VecType::Character : {
                SEXP charVec = PROTECT(Rf_duplicate(Rv));
                SEXP vectorPass = PROTECT(Rf_allocVector(STRSXP, m));

                SEXP res = ApplyFunction(charVec, vectorPass, mySample,
                                         myVec.get(), myReps, stdFun, myEnv,
                                         RFunVal, nthResFun, m, sampSize,
                                         IsNamed, IsGmp, n);
                UNPROTECT(2);
                return res;
            } case VecType::Complex : {
                SEXP vectorPass = PROTECT(Rf_allocVector(CPLXSXP, m));
                Rcomplex* ptr_vec = COMPLEX(vectorPass);

                Rcomplex* cmplxVec = COMPLEX(Rv);
                std::vector<Rcomplex> vCmplx(cmplxVec, cmplxVec + n);

                SEXP res = ApplyFunction(vCmplx, vectorPass, ptr_vec, mySample,
                                         myVec.get(), myReps, stdFun, myEnv,
                                         RFunVal, nthResFun, m, sampSize,
                                         IsNamed, IsGmp, n);
                UNPROTECT(1);
                return res;
            } case VecType::Raw : {
                SEXP vectorPass = PROTECT(Rf_allocVector(RAWSXP, m));
                Rbyte* ptr_vec = RAW(vectorPass);

                Rbyte* rawVec = RAW(Rv);
                std::vector<Rbyte> vByte(rawVec, rawVec + n);

                SEXP res = ApplyFunction(vByte, vectorPass, ptr_vec, mySample,
                                         myVec.get(), myReps, stdFun, myEnv,
                                         RFunVal, nthResFun, m, sampSize,
                                         IsNamed, IsGmp, n);
                UNPROTECT(1);
                return res;
            } case VecType::Logical : {
                SEXP vectorPass = PROTECT(Rf_allocVector(LGLSXP, m));
                int* ptr_vec = LOGICAL(vectorPass);

                SEXP res = ApplyFunction(vInt, vectorPass, ptr_vec, mySample,
                                         myVec.get(), myReps, stdFun, myEnv,
                                         RFunVal, nthResFun, m, sampSize,
                                         IsNamed, IsGmp, n);
                UNPROTECT(1);
                return res;
            } case VecType::Integer : {
                SEXP vectorPass = PROTECT(Rf_allocVector(INTSXP, m));
                int* ptr_vec = INTEGER(vectorPass);

                SEXP res = ApplyFunction(vInt, vectorPass, ptr_vec, mySample,
                                         myVec.get(), myReps, stdFun, myEnv,
                                         RFunVal, nthResFun, m, sampSize,
                                         IsNamed, IsGmp, n);
                UNPROTECT(1);
                return res;
            } default : {
                SEXP vectorPass = PROTECT(Rf_allocVector(REALSXP, m));
                double* ptr_vec = REAL(vectorPass);
                SEXP res = ApplyFunction(vNum, vectorPass, ptr_vec, mySample,
                                         myVec.get(), myReps, stdFun, myEnv,
                                         RFunVal, nthResFun, m, sampSize,
                                         IsNamed, IsGmp, n);
                UNPROTECT(1);
                return res;
            }
        }
    }

    switch (myType) {
        case VecType::Character : {
            SEXP charVec = PROTECT(Rf_duplicate(Rv));
            SEXP res = PROTECT(Rf_allocMatrix(STRSXP, sampSize, m));
            
            SampleResults(res, charVec, mySample, myVec.get(), myReps,
                          nthResFun, m, sampSize, n, IsGmp, IsNamed);
            
            UNPROTECT(2);
            return res;
        } case VecType::Complex : {
            std::vector<Rcomplex> stlCmplxVec(n);
            Rcomplex* vecCmplx = COMPLEX(Rv);

            for (int i = 0; i < n; ++i) {
                stlCmplxVec[i] = vecCmplx[i];
            }

            SEXP res = PROTECT(Rf_allocMatrix(CPLXSXP, sampSize, m));
            Rcomplex* matCmplx = COMPLEX(res);

            SampNoThrdSafe(matCmplx, res, stlCmplxVec, mySample,
                           myVec.get(), myReps, nthResFun, m, sampSize,
                           n, IsGmp, IsNamed);
            UNPROTECT(1);
            return res;
        } case VecType::Raw : {
            std::vector<Rbyte> stlRawVec(n);
            Rbyte* rawVec = RAW(Rv);

            for (int i = 0; i < n; ++i) {
                stlRawVec[i] = rawVec[i];
            }

            SEXP res = PROTECT(Rf_allocMatrix(RAWSXP, sampSize, m));
            Rbyte* rawMat = RAW(res);

            SampNoThrdSafe(rawMat, res, stlRawVec, mySample,
                           myVec.get(), myReps, nthResFun, m, sampSize,
                           n, IsGmp, IsNamed);
            UNPROTECT(1);
            return res;
        } case VecType::Logical : {
            vInt.assign(n, 0);
            int* vecBool = LOGICAL(Rv);

            for (int i = 0; i < n; ++i) {
                vInt[i] = vecBool[i];
            }

            SEXP res = PROTECT(Rf_allocMatrix(LGLSXP, sampSize, m));
            int* matBool = LOGICAL(res);

            SampNoThrdSafe(matBool, res, vInt, mySample, myVec.get(), myReps,
                           nthResFun, m, sampSize, n, IsGmp, IsNamed);
            UNPROTECT(1);
            return res;
        } case VecType::Integer : {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, sampSize, m));
            int* matInt = INTEGER(res);

            ThreadSafeSample(matInt, res, vInt, mySample, myVec.get(),
                             myReps, nthResFun, m, sampSize, nThreads,
                             Parallel, IsNamed, IsGmp, n);

            if (Rf_isFactor(Rv)) {
                SetFactorClass(res, Rv);
            }

            UNPROTECT(1);
            return res;
        } default : {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, sampSize, m));
            double* matNum = REAL(res);

            ThreadSafeSample(matNum, res, vNum, mySample, myVec.get(),
                             myReps, nthResFun, m, sampSize, nThreads,
                             Parallel, IsNamed, IsGmp, n);

            UNPROTECT(1);
            return res;
        }
    }
}
