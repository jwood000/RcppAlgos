#include "SetUpUtils.h"
#include "NthResult.h"
#include "RMatrix.h"
#include <thread>

template <typename T>
void SampleResults(T* sampleMatrix, const std::vector<T> &v,
                   const std::vector<double> &mySample,
                   mpz_t *const myBigSamp, const std::vector<int> &myReps,
                   nthResultPtr nthResFun, int m, int sampSize,
                   int lenV, bool IsGmp) {

    if (IsGmp) {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, 0.0,
                                                 myBigSamp[i], myReps);

            for (int j = 0; j < m; ++j) {
                sampleMatrix[i + sampSize * j] = v[z[j]];
            }
        }
    } else {
        mpz_t mpzDefault;
        mpz_init(mpzDefault);

        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[i],
                                                 mpzDefault, myReps);

            for (int j = 0; j < m; ++j) {
                sampleMatrix[i + sampSize * j] = v[z[j]];
            }
        }

        mpz_clear(mpzDefault);
    }
}

template <typename T>
void SampleResults(RcppParallel::RMatrix<T> &sampleMatrix,
                   const std::vector<T> &v,
                   const std::vector<double> &mySample,
                   mpz_t *const myBigSamp, const std::vector<int> &myReps,
                   nthResultPtr nthResFun, int m, int strtIdx, int endIdx,
                   int lenV, bool IsGmp) {

    if (IsGmp) {
        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, 0.0,
                                                 myBigSamp[i], myReps);

            for (int j = 0; j < m; ++j) {
                sampleMatrix(i, j) = v[z[j]];
            }
        }
    } else {
        mpz_t mpzDefault;
        mpz_init(mpzDefault);

        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[i],
                                                 mpzDefault, myReps);

            for (int j = 0; j < m; ++j) {
                sampleMatrix(i, j) = v[z[j]];
            }
        }

        mpz_clear(mpzDefault);
    }
}

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

void SampleResults(SEXP sampleMatrix, SEXP v,
                   const std::vector<double> &mySample,
                   mpz_t *const myBigSamp, const std::vector<int> &myReps,
                   nthResultPtr nthResFun, int m, int sampSize,
                   int lenV, bool IsGmp, bool IsNamed) {

    if (IsGmp) {
        for (int count = 0; count < sampSize; ++count) {
            const std::vector<int> z = nthResFun(lenV, m, 0.0,
                                                 myBigSamp[count], myReps);


            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(sampleMatrix, count + j * sampSize,
                               STRING_ELT(v, z[j]));
            }
        }
    } else {
        mpz_t mpzDefault;
        mpz_init(mpzDefault);

        for (int count = 0; count < sampSize; ++count) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[count],
                                                 mpzDefault, myReps);


            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(sampleMatrix, count + j * sampSize,
                               STRING_ELT(v, z[j]));
            }
        }

        mpz_clear(mpzDefault);
    }

    if (IsNamed) {
        SetSampleNames(sampleMatrix, IsGmp, sampSize, mySample, myBigSamp);
    }
}

SEXP SampCombPermMain(SEXP Rv, const std::vector<int> &vInt,
                      const std::vector<double> &vNum,
                      const std::vector<double> &mySample,
                      mpz_t *const myBigSamp,
                      const std::vector<int> &myReps,
                      nthResultPtr nthResFun, VecType myType, int n,
                      int m, int sampSize, int nThreads, bool IsNamed,
                      bool IsGmp, bool Parallel) {

    switch (myType) {
        case VecType::Character : {
            SEXP charVec = PROTECT(Rf_duplicate(Rv));
            SEXP res = PROTECT(Rf_allocMatrix(STRSXP, sampSize, m));

            SampleResults(res, charVec, mySample, myBigSamp, myReps,
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
                           myBigSamp, myReps, nthResFun, m, sampSize,
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
                           myBigSamp, myReps, nthResFun, m, sampSize,
                           n, IsGmp, IsNamed);
            UNPROTECT(1);
            return res;
        } case VecType::Logical : {
            std::vector<int> vBool(n, 0);
            int* vecBool = LOGICAL(Rv);

            for (int i = 0; i < n; ++i) {
                vBool[i] = vecBool[i];
            }

            SEXP res = PROTECT(Rf_allocMatrix(LGLSXP, sampSize, m));
            int* matBool = LOGICAL(res);

            SampNoThrdSafe(matBool, res, vBool, mySample, myBigSamp, myReps,
                           nthResFun, m, sampSize, n, IsGmp, IsNamed);
            UNPROTECT(1);
            return res;
        } case VecType::Integer : {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, sampSize, m));
            int* matInt = INTEGER(res);

            ThreadSafeSample(matInt, res, vInt, mySample, myBigSamp,
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

            ThreadSafeSample(matNum, res, vNum, mySample, myBigSamp,
                             myReps, nthResFun, m, sampSize, nThreads,
                             Parallel, IsNamed, IsGmp, n);

            UNPROTECT(1);
            return res;
        }
    }
}
