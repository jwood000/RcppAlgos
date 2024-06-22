#include "SetUpUtils.h"
#include "NthResult.h"
#include "RMatrix.h"
#include <thread>

template <typename T>
void SampleResults(T* sampleMatrix, const std::vector<T> &v,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   const std::vector<int> &myReps,
                   nthResultPtr nthResFun, std::size_t m,
                   std::size_t sampSize, int lenV, bool IsGmp) {

    if (IsGmp) {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, 0.0,
                                                 myBigSamp[i], myReps);

            for (std::size_t j = 0; j < m; ++j) {
                sampleMatrix[i + sampSize * j] = v[z[j]];
            }
        }
    } else {
        mpz_class mpzDefault;

        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[i],
                                                 mpzDefault, myReps);

            for (std::size_t j = 0; j < m; ++j) {
                sampleMatrix[i + sampSize * j] = v[z[j]];
            }
        }
    }
}

template <typename T>
void SampleResults(RcppParallel::RMatrix<T> &sampleMatrix,
                   const std::vector<T> &v,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   const std::vector<int> &myReps, nthResultPtr nthResFun,
                   int m, int strtIdx, int endIdx, int lenV, bool IsGmp) {

    if (IsGmp) {
        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, 0.0,
                                                 myBigSamp[i], myReps);

            for (int j = 0; j < m; ++j) {
                sampleMatrix(i, j) = v[z[j]];
            }
        }
    } else {
        mpz_class mpzDefault;

        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[i],
                                                 mpzDefault, myReps);

            for (int j = 0; j < m; ++j) {
                sampleMatrix(i, j) = v[z[j]];
            }
        }
    }
}

template <typename T>
void SampNoThrdSafe(T* sampleMatrix, SEXP res, const std::vector<T> &v,
                    const std::vector<double> &mySample,
                    const std::vector<mpz_class> &myBigSamp,
                    const std::vector<int> &myReps,
                    nthResultPtr nthResFun, int m, int sampSize,
                    int lenV, bool IsGmp, bool IsNamed) {

    SampleResults(sampleMatrix, v, mySample, myBigSamp, myReps,
                  nthResFun, m, sampSize, lenV, IsGmp);
    SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp, IsNamed);
}

template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T> &sampleMatrix,
                  const std::vector<T> &v,
                  const std::vector<double> &mySample,
                  const std::vector<mpz_class> &myBigSamp,
                  const std::vector<int> &myReps,
                  nthResultPtr nthResFun, int m, int strtIdx, int endIdx,
                  int lenV, bool IsGmp) {

    SampleResults(sampleMatrix, v, mySample, myBigSamp, myReps,
                  nthResFun, m, strtIdx, endIdx, lenV, IsGmp);
}

template <typename T>
void ThreadSafeSample(T* mat, SEXP res, const std::vector<T> &v,
                      const std::vector<double> &mySample,
                      const std::vector<mpz_class> &myBigSamp,
                      const std::vector<int> &myReps,
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
                                 std::cref(mySample), std::cref(myBigSamp),
                                 std::cref(myReps), nthResFun,
                                 m, step, nextStep, lenV, IsGmp);
        }

        threads.emplace_back(
            std::cref(ParallelGlue<T>), std::ref(parMat), std::cref(v),
            std::cref(mySample),std::cref(myBigSamp), std::cref(myReps),
            nthResFun, m, step, sampSize, lenV, IsGmp
        );

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        SampleResults(mat, v, mySample, myBigSamp, myReps,
                      nthResFun, m, sampSize, lenV, IsGmp);
    }

    SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp, IsNamed);
}

void SampleResults(SEXP sampleMatrix, SEXP v,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   const std::vector<int> &myReps,
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
        mpz_class mpzDefault;

        for (int count = 0; count < sampSize; ++count) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[count],
                                                 mpzDefault, myReps);


            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(sampleMatrix, count + j * sampSize,
                               STRING_ELT(v, z[j]));
            }
        }
    }

    SetSampleNames(sampleMatrix, IsGmp, sampSize,
                   mySample, myBigSamp, IsNamed);
}

SEXP SampCombPermMain(SEXP Rv, const std::vector<int> &vInt,
                      const std::vector<double> &vNum,
                      const std::vector<double> &mySample,
                      const std::vector<mpz_class> &myBigSamp,
                      const std::vector<int> &myReps,
                      nthResultPtr nthResFun, VecType myType, int n,
                      int m, int sampSize, int nThreads, bool IsNamed,
                      bool IsGmp, bool Parallel) {

    switch (myType) {
        case VecType::Character : {
            cpp11::sexp charVec = Rf_duplicate(Rv);
            cpp11::sexp res = Rf_allocMatrix(STRSXP, sampSize, m);

            SampleResults(res, charVec, mySample, myBigSamp, myReps,
                          nthResFun, m, sampSize, n, IsGmp, IsNamed);

            return res;
        } case VecType::Complex : {
            std::vector<Rcomplex> stlCmplxVec =
                CppConvert::GetVec<Rcomplex>(Rv);
            cpp11::sexp res = Rf_allocMatrix(CPLXSXP, sampSize, m);
            Rcomplex* matCmplx = COMPLEX(res);

            SampNoThrdSafe(matCmplx, res, stlCmplxVec, mySample,
                           myBigSamp, myReps, nthResFun, m, sampSize,
                           n, IsGmp, IsNamed);

            return res;
        } case VecType::Raw : {
            std::vector<Rbyte> stlRawVec = CppConvert::GetVec<Rbyte>(Rv);
            cpp11::sexp res = Rf_allocMatrix(RAWSXP, sampSize, m);
            Rbyte* rawMat = RAW(res);

            SampNoThrdSafe(rawMat, res, stlRawVec, mySample,
                           myBigSamp, myReps, nthResFun, m, sampSize,
                           n, IsGmp, IsNamed);

            return res;
        } case VecType::Logical : {
            std::vector<int> vBool = CppConvert::GetVec<int>(Rv);
            cpp11::sexp res = Rf_allocMatrix(LGLSXP, sampSize, m);
            int* matBool = LOGICAL(res);

            SampNoThrdSafe(matBool, res, vBool, mySample, myBigSamp, myReps,
                           nthResFun, m, sampSize, n, IsGmp, IsNamed);

            return res;
        } case VecType::Integer : {
            cpp11::sexp res = Rf_allocMatrix(INTSXP, sampSize, m);
            int* matInt = INTEGER(res);

            ThreadSafeSample(matInt, res, vInt, mySample, myBigSamp,
                             myReps, nthResFun, m, sampSize, nThreads,
                             Parallel, IsNamed, IsGmp, n);

            if (Rf_isFactor(Rv)) {
                SetFactorClass(res, Rv);
            }

            return res;
        } default : {
            cpp11::sexp res = Rf_allocMatrix(REALSXP, sampSize, m);
            double* matNum = REAL(res);

            ThreadSafeSample(matNum, res, vNum, mySample, myBigSamp,
                             myReps, nthResFun, m, sampSize, nThreads,
                             Parallel, IsNamed, IsGmp, n);

            return res;
        }
    }
}
