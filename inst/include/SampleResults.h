#ifndef SAMPLE_RESULTS_H
#define SAMPLE_RESULTS_H

#include "SetUpUtils.h"
#include "NthResult.h"
#include "RMatrix.h"

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

#endif
