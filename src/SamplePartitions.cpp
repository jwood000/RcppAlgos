#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/NthPartition.h"
#include "RMatrix.h"
#include <thread>

template <typename T>
void SampleResults(T* sampleMatrix, const std::vector<T> &v,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   const std::vector<int> &myReps,
                   nthPartsPtr nthPartFun, std::size_t m,
                   std::size_t sampSize, int tar, int strtLen,
                   int cap, bool IsGmp) {

    if (IsGmp) {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthPartFun(tar, m, cap, strtLen,
                                                  0.0, myBigSamp[i]);
            for (std::size_t j = 0; j < m; ++j) {
                sampleMatrix[i + sampSize * j] = v[z[j]];
            }
        }
    } else {
        mpz_class mpzDefault;

        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthPartFun(tar, m, cap, strtLen,
                                                  mySample[i], mpzDefault);

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
                   const std::vector<int> &myReps,
                   nthPartsPtr nthPartFun, int m, int strtIdx, int endIdx,
                   int tar, int strtLen, int cap, bool IsGmp) {

    if (IsGmp) {
        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthPartFun(tar, m, cap, strtLen,
                                                  0.0, myBigSamp[i]);
            for (int j = 0; j < m; ++j) {
                sampleMatrix(i, j) = v[z[j]];
            }
        }
    } else {
        mpz_class mpzDefault;

        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthPartFun(tar, m, cap, strtLen,
                                                  mySample[i], mpzDefault);
            for (int j = 0; j < m; ++j) {
                sampleMatrix(i, j) = v[z[j]];
            }
        }
    }
}

template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T> &sampleMatrix,
                  const std::vector<T> &v,
                  const std::vector<double> &mySample,
                  const std::vector<mpz_class> &myBigSamp,
                  const std::vector<int> &myReps,
                  nthPartsPtr nthPartFun, int m, int strtIdx, int endIdx,
                  int tar, int strtLen, int cap, bool IsGmp) {

    SampleResults(sampleMatrix, v, mySample, myBigSamp, myReps,
                  nthPartFun, m, strtIdx, endIdx, tar, strtLen, cap, IsGmp);
}

template <typename T>
void ThreadSafeSample(T* mat, SEXP res, const std::vector<T> &v,
                      const std::vector<double> &mySample,
                      const std::vector<mpz_class> &myBigSamp,
                      const std::vector<int> &myReps,
                      nthPartsPtr nthPartFun, int m, int sampSize,
                      int nThreads, bool Parallel, bool IsNamed,
                      int tar, int strtLen, int cap, bool IsGmp) {

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
                                 std::cref(myReps), nthPartFun, m,
                                 step, nextStep, tar, strtLen, cap, IsGmp);
        }

        threads.emplace_back(
            std::cref(ParallelGlue<T>), std::ref(parMat), std::cref(v),
            std::cref(mySample), std::cref(myBigSamp), std::cref(myReps),
            nthPartFun, m, step, sampSize, tar, strtLen, cap, IsGmp
        );

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        SampleResults(mat, v, mySample, myBigSamp, myReps,
                      nthPartFun, m, sampSize, tar, strtLen, cap, IsGmp);
    }

    SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp, IsNamed);
}

[[cpp11::register]]
SEXP SamplePartitions(
    SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP RindexVec,
    SEXP RmySeed, SEXP RNumSamp, SEXP baseSample, SEXP RNumThreads,
    SEXP RmaxThreads, SEXP RNamed, SEXP RcompFun, SEXP Rtarget,
    SEXP myEnv, SEXP RIsComposition, SEXP RIsWeak
) {

    int n = 0;
    int m = 0;
    int nThreads = 1;
    int maxThreads = 1;

    bool IsMult = false;
    VecType myType = VecType::Integer;
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");

    std::vector<double> vNum;
    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;

    bool Parallel = false; // This will be set in SetThreads below. For the
          // partition and composition functions we don't have a Parallel
          // argument. The goal is to eventually phase out this argument.
    bool IsRep    = CppConvert::convertFlag(RisRep, "repetition");
    bool IsNamed  = CppConvert::convertFlag(RNamed, "namedSample");

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum, Rv,
              RFreqs, Rm, n, m, IsMult, IsRep, true);

    const std::string mainFun = "sum";

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> targetIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compVec;
    std::vector<double> targetVals;

    ConstraintType ctype = ConstraintType::NoConstraint;
    PartDesign part;

    part.isRep   = IsRep;
    part.isMult  = IsMult;
    part.mIsNull = Rf_isNull(Rm);
    part.isWeak  = CppConvert::convertFlag(RIsWeak, "weak");
    part.isComp  = CppConvert::convertFlag(RIsComposition,
                                             "IsComposition");
    part.isComb = !part.isComp;

    cpp11::sexp Rlow = R_NilValue;
    ConstraintSetup(vNum, myReps, targetVals, vInt, targetIntVals, funDbl,
                    part, ctype, n, m, compVec, mainFun, mainFun, myType,
                    Rtarget, RcompFun, R_NilValue, Rlow);

    if (part.ptype == PartitionType::Multiset ||
        part.ptype == PartitionType::CoarseGrained ||
        part.ptype == PartitionType::NotPartition) {

        cpp11::stop("Partition sampling not available for this case.");
    }

    int sampSize;
    std::vector<double> mySample;

    // This can occur if we are dealing with capped cases where calculating
    // the number of partitions could take a long time. When this occurs with
    // partitionsGeneral, it is faster to generate partitions and push them
    // to a vector until the next partitions algorithm exhaust, then we can
    // convert this to an R matrix (instead of preallocating a matrix).
    //
    // When we are dealing with sampling, we have to know the total number
    // of partitions, thus the following:

    if (part.numUnknown) PartitionsCount(myReps, part, n, true);
    const bool SampleGmp = (part.count > SampleLimit);

    if (SampleGmp && !part.isGmp) {
        part.isGmp    = true;
        part.bigCount = part.count;
    }

    SetRandomSample(RindexVec, RNumSamp, sampSize, SampleGmp,
                    part.count, mySample, baseSample, myEnv);

    const int bigSampSize = (SampleGmp) ? sampSize : 1;
    std::vector<mpz_class> myVec(bigSampSize);

    SetRandomSampleMpz(RindexVec, RmySeed, sampSize,
                       SampleGmp, part.bigCount, myVec);

    const int limit = 2;
    SetThreads(Parallel, maxThreads, sampSize,
               myType, nThreads, RNumThreads, limit);

    const int cap     = n - static_cast<int>(part.includeZero);
    const int strtLen = std::count_if(part.startZ.cbegin(),
                                      part.startZ.cend(),
                                      [](int i){return i > 0;});

    if (myType == VecType::Integer) {
        cpp11::sexp res = Rf_allocMatrix(INTSXP, sampSize, part.width);
        int* matInt = INTEGER(res);

        if (part.width == 1) {
            matInt[0] = Rf_asInteger(Rtarget);
            SetSampleNames(res, false, sampSize,
                           mySample, myVec, IsNamed);
        } else {
            const nthPartsPtr nthPartFun = GetNthPartsFunc(
                part.ptype, part.isGmp, part.isComp
            );

            ThreadSafeSample(matInt, res, vInt, mySample, myVec,
                             myReps, nthPartFun, part.width, sampSize,
                             nThreads, Parallel, IsNamed, part.mapTar,
                             strtLen, cap, part.isGmp);
        }

        return res;
    } else {
        cpp11::sexp res = Rf_allocMatrix(REALSXP, sampSize, part.width);
        double* matNum = REAL(res);

        if (part.width == 1) {
            matNum[0] = Rf_asReal(Rtarget);
            SetSampleNames(res, false, sampSize,
                           mySample, myVec, IsNamed);
        } else {
            const nthPartsPtr nthPartFun = GetNthPartsFunc(
                part.ptype, part.isGmp, part.isComp
            );

            ThreadSafeSample(matNum, res, vNum, mySample, myVec,
                             myReps, nthPartFun, part.width, sampSize,
                             nThreads, Parallel, IsNamed, part.mapTar,
                             strtLen, cap, part.isGmp);
        }

        return res;
    }
}

template void ThreadSafeSample(int*, SEXP, const std::vector<int>&,
                               const std::vector<double>&,
                               const std::vector<mpz_class>&,
                               const std::vector<int>&, nthPartsPtr, int,
                               int, int, bool, bool, int, int, int, bool);

template void ThreadSafeSample(double*, SEXP, const std::vector<double>&,
                               const std::vector<double>&,
                               const std::vector<mpz_class>&,
                               const std::vector<int>&, nthPartsPtr, int,
                               int, int, bool, bool, int, int, int, bool);
