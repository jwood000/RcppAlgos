#include "Sample/SampCombPermStd.h"
#include "Sample/SampleCombPerm.h"
#include "Sample/SampleApply.h"
#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "SetUpUtils.h"

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
    bool IsNamed = CleanConvert::convertFlag(RNamed, "namedSample");

    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;
    std::vector<double> vNum;

    bool Parallel = CleanConvert::convertFlag(Rparallel, "Parallel");
    bool IsRep = CleanConvert::convertFlag(RisRep, "repetition");
    const bool IsComb = CleanConvert::convertFlag(RIsComb, "IsComb");
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

        return SampleCombPermApply(Rv, vInt, vNum, mySample, myVec.get(),
                                   myReps, stdFun, myEnv, RFunVal, nthResFun,
                                   myType, n, m, sampSize, IsNamed, IsGmp);
    }

    return SampCombPermMain(Rv, vInt, vNum, mySample, myVec.get(),
                            myReps, nthResFun, myType, n, m, sampSize,
                            nThreads, IsNamed, IsGmp, Parallel);
}
