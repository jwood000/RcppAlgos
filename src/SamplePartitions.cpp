#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/NthPartition.h"
#include "SamplePartitions.h"
#include "Cpp14MakeUnique.h"
#include "RMatrix.h"
#include <thread>

template <typename T>
void SampleResults(T* sampleMatrix, const std::vector<T> &v,
                   const std::vector<double> &mySample,
                   mpz_t *const myBigSamp, const std::vector<int> &myReps,
                   nthPartsPtr nthPartFun, int m, int sampSize,
                   int tar, int k, bool IsGmp) {

    if (IsGmp) {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthPartFun(tar, m, k, 0.0,
                                                  myBigSamp[i]);

            for (int j = 0; j < m; ++j) {
                sampleMatrix[i + sampSize * j] = v[z[j]];
            }
        }
    } else {
        mpz_t mpzDefault;
        mpz_init(mpzDefault);

        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthPartFun(tar, m, k,
                                                  mySample[i],
                                                  mpzDefault);

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
                   nthPartsPtr nthPartFun, int m, int strtIdx, int endIdx,
                   int tar, int k, bool IsGmp) {

    if (IsGmp) {
        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthPartFun(tar, m, k, 0.0,
                                                  myBigSamp[i]);
            for (int j = 0; j < m; ++j) {
                sampleMatrix(i, j) = v[z[j]];
            }
        }
    } else {
        mpz_t mpzDefault;
        mpz_init(mpzDefault);

        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthPartFun(tar, m, k,
                                                  mySample[i],
                                                  mpzDefault);
            for (int j = 0; j < m; ++j) {
                sampleMatrix(i, j) = v[z[j]];
            }
        }

        mpz_clear(mpzDefault);
    }
}

template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T> &sampleMatrix,
                  const std::vector<T> &v,
                  const std::vector<double> &mySample,
                  mpz_t *const myBigSamp, const std::vector<int> &myReps,
                  nthPartsPtr nthPartFun, int m, int strtIdx, int endIdx,
                  int tar, int k, bool IsGmp) {

    SampleResults(sampleMatrix, v, mySample, myBigSamp, myReps,
                  nthPartFun, m, strtIdx, endIdx, tar, k, IsGmp);
}

template <typename T>
void ThreadSafeSample(T* mat, SEXP res, const std::vector<T> &v,
                      const std::vector<double> &mySample,
                      mpz_t *const myBigSamp, const std::vector<int> &myReps,
                      nthPartsPtr nthPartFun, int m, int sampSize,
                      int nThreads, bool Parallel, bool IsNamed,
                      int tar, int k, bool IsGmp) {

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
                                 std::cref(myReps), nthPartFun, m,
                                 step, nextStep, tar, k, IsGmp);
        }

        threads.emplace_back(std::cref(ParallelGlue<T>), std::ref(parMat),
                             std::cref(v), std::cref(mySample), myBigSamp,
                             std::cref(myReps), nthPartFun, m, step,
                             sampSize, tar, k, IsGmp);

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        SampleResults(mat, v, mySample, myBigSamp, myReps,
                      nthPartFun, m, sampSize, tar, k, IsGmp);
    }

    if (IsNamed) {
        SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp);
    }
}

SEXP SamplePartitions(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                      SEXP RindexVec, SEXP RmySeed, SEXP RNumSamp,
                      SEXP baseSample, SEXP Rparallel, SEXP RNumThreads,
                      SEXP RmaxThreads, SEXP RNamed, SEXP RcompFun,
                      SEXP Rtarget, SEXP Rtolerance, SEXP myEnv) {

    int n = 0;
    int m = 0;
    int nThreads = 1;
    int maxThreads = 1;

    bool IsMult = false;
    VecType myType = VecType::Integer;
    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

    std::vector<double> vNum;
    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;

    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRep    = CleanConvert::convertLogical(RisRep, "repetition");
    bool IsNamed  = CleanConvert::convertLogical(RNamed, "namedSample");

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum, Rv,
              RFreqs, Rm, n, m, IsMult, IsRep, true);

    const std::string mainFun = "sum";

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> targetIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compFunVec;
    std::vector<double> targetVals;

    ConstraintType ctype = ConstraintType::NoConstraint;
    PartDesign part;

    part.isRep = IsRep;
    part.isMult = IsMult;
    part.mIsNull = Rf_isNull(Rm);
    SEXP Rlow = R_NilValue;

    ConstraintSetup(vNum, myReps, targetVals, targetIntVals, funDbl,
                    part, ctype, n, m, compFunVec, mainFun, myType,
                    Rtarget, RcompFun, Rtolerance, Rlow, true, false);

    if (part.ptype == PartitionType::Multiset ||
        part.ptype == PartitionType::CoarseGrained ||
        part.ptype == PartitionType::NotPartition) {

        Rf_error("Partition sampling not available for this case.");
    }

    int sampSize;
    std::vector<double> mySample;
    SetRandomSample(RindexVec, RNumSamp, sampSize, part.isGmp,
                    part.count, mySample, baseSample, myEnv);

    const int bigSampSize = (part.isGmp) ? sampSize : 1;
    auto myVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);

    for (int i = 0; i < bigSampSize; ++i) {
        mpz_init(myVec[i]);
    }

    SetRandomSampleMpz(RindexVec, RmySeed, sampSize,
                       part.isGmp, part.bigCount, myVec.get());

    const int limit = 2;
    SetThreads(Parallel, maxThreads, sampSize,
               myType, nThreads, RNumThreads, limit);

    const int strtLen = std::count_if(part.startZ.cbegin(),
                                      part.startZ.cend(),
                                      [](int i){return i > 0;});

    const int k = (part.ptype == PartitionType::DstctSpecial ||
                   part.ptype == PartitionType::DstctStdAll) ? strtLen : n;

    const nthPartsPtr nthPartFun = GetNthPartsFunc(part.ptype, part.isGmp);

    if (myType == VecType::Integer) {
        SEXP res = PROTECT(Rf_allocMatrix(INTSXP, sampSize, m));
        int* matInt = INTEGER(res);

        ThreadSafeSample(matInt, res, vInt, mySample, myVec.get(),
                         myReps, nthPartFun, part.width, sampSize, nThreads,
                         Parallel, IsNamed, part.mapTar, k, part.isGmp);

        UNPROTECT(1);
        return res;
    } else {
        SEXP res = PROTECT(Rf_allocMatrix(REALSXP, sampSize, m));
        double* matNum = REAL(res);

        ThreadSafeSample(matNum, res, vNum, mySample, myVec.get(),
                         myReps, nthPartFun, part.width, sampSize, nThreads,
                         Parallel, IsNamed, part.mapTar, k, part.isGmp);

        UNPROTECT(1);
        return res;
    }
}
