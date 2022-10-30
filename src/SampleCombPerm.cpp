#include "Sample/SampCombPermStd.h"
#include "Sample/SampleApply.h"
#include "ComputedCount.h"
#include "SetUpUtils.h"

[[cpp11::register]]
SEXP SampleCombPerm(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                    SEXP RindexVec, SEXP RIsComb, SEXP RmySeed,
                    SEXP RNumSamp, SEXP baseSample, SEXP stdFun, SEXP myEnv,
                    SEXP Rparallel, SEXP RNumThreads, SEXP RmaxThreads,
                    SEXP RNamed, SEXP RFunVal) {

    int n = 0;
    int m = 0;
    int nThreads = 1;
    int maxThreads = 1;
    bool applyFun = false;

    if (!Rf_isNull(stdFun) && !Rf_isFactor(Rv)) {
        if (!Rf_isFunction(stdFun)) {
            cpp11::stop("FUN must be a function!");
        }

        applyFun = true;
    }

    VecType myType = VecType::Integer;
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");
    bool IsNamed = CppConvert::convertFlag(RNamed, "namedSample");

    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;
    std::vector<double> vNum;

    bool Parallel = CppConvert::convertFlag(Rparallel, "Parallel");
    bool IsRep = CppConvert::convertFlag(RisRep, "repetition");
    const bool IsComb = CppConvert::convertFlag(RIsComb, "IsComb");
    bool IsMult = false;

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum,
              Rv, RFreqs, Rm, n, m, IsMult, IsRep);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp = (computedRows > SampleLimit);

    mpz_class computedRowsMpz;

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult, IsComb,
                          IsRep, n, m, Rm, freqs, myReps);
    }

    int sampSize;
    std::vector<double> mySample;
    SetRandomSample(RindexVec, RNumSamp, sampSize, IsGmp,
                    computedRows, mySample, baseSample, myEnv);

    const int bigSampSize = IsGmp ? sampSize : 1;
    std::vector<mpz_class> myVec(bigSampSize);

    SetRandomSampleMpz(RindexVec, RmySeed, sampSize,
                       IsGmp, computedRowsMpz, myVec);

    const int limit = 2;
    SetThreads(Parallel, maxThreads, sampSize,
               myType, nThreads, RNumThreads, limit);

    const nthResultPtr nthResFun = GetNthResultFunc(IsComb, IsMult,
                                                    IsRep, IsGmp);

    if (applyFun) {
        return SampleCombPermApply(Rv, vInt, vNum, mySample, myVec,
                                   myReps, stdFun, myEnv, RFunVal, nthResFun,
                                   myType, n, m, sampSize, IsNamed, IsGmp);
    }

    return SampCombPermMain(Rv, vInt, vNum, mySample, myVec,
                            myReps, nthResFun, myType, n, m, sampSize,
                            nThreads, IsNamed, IsGmp, Parallel);
}
