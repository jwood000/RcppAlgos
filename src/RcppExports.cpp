// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// CheckReturn
int CheckReturn(SEXP Rv, SEXP f1, SEXP f2, SEXP Rtarget, bool IsFactor, SEXP RKeepRes, SEXP stdFun);
RcppExport SEXP _RcppAlgos_CheckReturn(SEXP RvSEXP, SEXP f1SEXP, SEXP f2SEXP, SEXP RtargetSEXP, SEXP IsFactorSEXP, SEXP RKeepResSEXP, SEXP stdFunSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rtarget(RtargetSEXP);
    Rcpp::traits::input_parameter< bool >::type IsFactor(IsFactorSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RKeepRes(RKeepResSEXP);
    Rcpp::traits::input_parameter< SEXP >::type stdFun(stdFunSEXP);
    rcpp_result_gen = Rcpp::wrap(CheckReturn(Rv, f1, f2, Rtarget, IsFactor, RKeepRes, stdFun));
    return rcpp_result_gen;
END_RCPP
}
// CombinatoricsStndrd
SEXP CombinatoricsStndrd(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow, SEXP Rhigh, bool IsComb, bool IsFactor, SEXP Rparallel, SEXP RNumThreads, int maxThreads);
RcppExport SEXP _RcppAlgos_CombinatoricsStndrd(SEXP RvSEXP, SEXP RmSEXP, SEXP RisRepSEXP, SEXP RFreqsSEXP, SEXP RlowSEXP, SEXP RhighSEXP, SEXP IsCombSEXP, SEXP IsFactorSEXP, SEXP RparallelSEXP, SEXP RNumThreadsSEXP, SEXP maxThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rm(RmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RisRep(RisRepSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RFreqs(RFreqsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rlow(RlowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rhigh(RhighSEXP);
    Rcpp::traits::input_parameter< bool >::type IsComb(IsCombSEXP);
    Rcpp::traits::input_parameter< bool >::type IsFactor(IsFactorSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rparallel(RparallelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumThreads(RNumThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxThreads(maxThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(CombinatoricsStndrd(Rv, Rm, RisRep, RFreqs, Rlow, Rhigh, IsComb, IsFactor, Rparallel, RNumThreads, maxThreads));
    return rcpp_result_gen;
END_RCPP
}
// CombinatoricsApply
SEXP CombinatoricsApply(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow, SEXP Rhigh, bool IsComb, SEXP stdFun, SEXP myEnv);
RcppExport SEXP _RcppAlgos_CombinatoricsApply(SEXP RvSEXP, SEXP RmSEXP, SEXP RisRepSEXP, SEXP RFreqsSEXP, SEXP RlowSEXP, SEXP RhighSEXP, SEXP IsCombSEXP, SEXP stdFunSEXP, SEXP myEnvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rm(RmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RisRep(RisRepSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RFreqs(RFreqsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rlow(RlowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rhigh(RhighSEXP);
    Rcpp::traits::input_parameter< bool >::type IsComb(IsCombSEXP);
    Rcpp::traits::input_parameter< SEXP >::type stdFun(stdFunSEXP);
    Rcpp::traits::input_parameter< SEXP >::type myEnv(myEnvSEXP);
    rcpp_result_gen = Rcpp::wrap(CombinatoricsApply(Rv, Rm, RisRep, RFreqs, Rlow, Rhigh, IsComb, stdFun, myEnv));
    return rcpp_result_gen;
END_RCPP
}
// CombinatoricsCount
SEXP CombinatoricsCount(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, bool IsComb);
RcppExport SEXP _RcppAlgos_CombinatoricsCount(SEXP RvSEXP, SEXP RmSEXP, SEXP RisRepSEXP, SEXP RFreqsSEXP, SEXP IsCombSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rm(RmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RisRep(RisRepSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RFreqs(RFreqsSEXP);
    Rcpp::traits::input_parameter< bool >::type IsComb(IsCombSEXP);
    rcpp_result_gen = Rcpp::wrap(CombinatoricsCount(Rv, Rm, RisRep, RFreqs, IsComb));
    return rcpp_result_gen;
END_RCPP
}
// ComboGroupsCountCpp
SEXP ComboGroupsCountCpp(SEXP Rv, SEXP RNumGroups);
RcppExport SEXP _RcppAlgos_ComboGroupsCountCpp(SEXP RvSEXP, SEXP RNumGroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumGroups(RNumGroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(ComboGroupsCountCpp(Rv, RNumGroups));
    return rcpp_result_gen;
END_RCPP
}
// ComboGroupsRcpp
SEXP ComboGroupsRcpp(SEXP Rv, SEXP RNumGroups, SEXP RRetType, SEXP Rlow, SEXP Rhigh, bool IsFactor, SEXP Rparallel, SEXP RNumThreads, int maxThreads, bool IsSample, SEXP RindexVec, SEXP RmySeed, SEXP RNumSamp, Rcpp::Function baseSample, SEXP RNamed);
RcppExport SEXP _RcppAlgos_ComboGroupsRcpp(SEXP RvSEXP, SEXP RNumGroupsSEXP, SEXP RRetTypeSEXP, SEXP RlowSEXP, SEXP RhighSEXP, SEXP IsFactorSEXP, SEXP RparallelSEXP, SEXP RNumThreadsSEXP, SEXP maxThreadsSEXP, SEXP IsSampleSEXP, SEXP RindexVecSEXP, SEXP RmySeedSEXP, SEXP RNumSampSEXP, SEXP baseSampleSEXP, SEXP RNamedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumGroups(RNumGroupsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RRetType(RRetTypeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rlow(RlowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rhigh(RhighSEXP);
    Rcpp::traits::input_parameter< bool >::type IsFactor(IsFactorSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rparallel(RparallelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumThreads(RNumThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxThreads(maxThreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type IsSample(IsSampleSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RindexVec(RindexVecSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RmySeed(RmySeedSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumSamp(RNumSampSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type baseSample(baseSampleSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNamed(RNamedSEXP);
    rcpp_result_gen = Rcpp::wrap(ComboGroupsRcpp(Rv, RNumGroups, RRetType, Rlow, Rhigh, IsFactor, Rparallel, RNumThreads, maxThreads, IsSample, RindexVec, RmySeed, RNumSamp, baseSample, RNamed));
    return rcpp_result_gen;
END_RCPP
}
// CombinatoricsCnstrt
SEXP CombinatoricsCnstrt(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow, SEXP Rhigh, SEXP f1, SEXP f2, SEXP Rtarget, bool IsComb, SEXP RKeepRes, SEXP Rparallel, SEXP RnThreads, int maxThreads, SEXP Rtolerance);
RcppExport SEXP _RcppAlgos_CombinatoricsCnstrt(SEXP RvSEXP, SEXP RmSEXP, SEXP RisRepSEXP, SEXP RFreqsSEXP, SEXP RlowSEXP, SEXP RhighSEXP, SEXP f1SEXP, SEXP f2SEXP, SEXP RtargetSEXP, SEXP IsCombSEXP, SEXP RKeepResSEXP, SEXP RparallelSEXP, SEXP RnThreadsSEXP, SEXP maxThreadsSEXP, SEXP RtoleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rm(RmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RisRep(RisRepSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RFreqs(RFreqsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rlow(RlowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rhigh(RhighSEXP);
    Rcpp::traits::input_parameter< SEXP >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rtarget(RtargetSEXP);
    Rcpp::traits::input_parameter< bool >::type IsComb(IsCombSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RKeepRes(RKeepResSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rparallel(RparallelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RnThreads(RnThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxThreads(maxThreadsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rtolerance(RtoleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(CombinatoricsCnstrt(Rv, Rm, RisRep, RFreqs, Rlow, Rhigh, f1, f2, Rtarget, IsComb, RKeepRes, Rparallel, RnThreads, maxThreads, Rtolerance));
    return rcpp_result_gen;
END_RCPP
}
// DivNumSieve
SEXP DivNumSieve(SEXP Rb1, SEXP Rb2, bool bDivSieve, SEXP RNamed, SEXP RNumThreads, int maxThreads);
RcppExport SEXP _RcppAlgos_DivNumSieve(SEXP Rb1SEXP, SEXP Rb2SEXP, SEXP bDivSieveSEXP, SEXP RNamedSEXP, SEXP RNumThreadsSEXP, SEXP maxThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rb1(Rb1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rb2(Rb2SEXP);
    Rcpp::traits::input_parameter< bool >::type bDivSieve(bDivSieveSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNamed(RNamedSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumThreads(RNumThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxThreads(maxThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(DivNumSieve(Rb1, Rb2, bDivSieve, RNamed, RNumThreads, maxThreads));
    return rcpp_result_gen;
END_RCPP
}
// GetClassVals
Rcpp::List GetClassVals(bool IsStdRet, SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, bool IsComb, bool IsFactor, SEXP stdFun);
RcppExport SEXP _RcppAlgos_GetClassVals(SEXP IsStdRetSEXP, SEXP RvSEXP, SEXP RmSEXP, SEXP RisRepSEXP, SEXP RFreqsSEXP, SEXP IsCombSEXP, SEXP IsFactorSEXP, SEXP stdFunSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type IsStdRet(IsStdRetSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rm(RmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RisRep(RisRepSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RFreqs(RFreqsSEXP);
    Rcpp::traits::input_parameter< bool >::type IsComb(IsCombSEXP);
    Rcpp::traits::input_parameter< bool >::type IsFactor(IsFactorSEXP);
    Rcpp::traits::input_parameter< SEXP >::type stdFun(stdFunSEXP);
    rcpp_result_gen = Rcpp::wrap(GetClassVals(IsStdRet, Rv, Rm, RisRep, RFreqs, IsComb, IsFactor, stdFun));
    return rcpp_result_gen;
END_RCPP
}
// MotleyContainer
SEXP MotleyContainer(SEXP Rb1, SEXP Rb2, bool isEuler, SEXP RNamed, SEXP RNumThreads, int maxThreads);
RcppExport SEXP _RcppAlgos_MotleyContainer(SEXP Rb1SEXP, SEXP Rb2SEXP, SEXP isEulerSEXP, SEXP RNamedSEXP, SEXP RNumThreadsSEXP, SEXP maxThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rb1(Rb1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rb2(Rb2SEXP);
    Rcpp::traits::input_parameter< bool >::type isEuler(isEulerSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNamed(RNamedSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumThreads(RNumThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxThreads(maxThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(MotleyContainer(Rb1, Rb2, isEuler, RNamed, RNumThreads, maxThreads));
    return rcpp_result_gen;
END_RCPP
}
// cpp11GetNumThreads
int cpp11GetNumThreads();
RcppExport SEXP _RcppAlgos_cpp11GetNumThreads() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(cpp11GetNumThreads());
    return rcpp_result_gen;
END_RCPP
}
// PollardRhoContainer
SEXP PollardRhoContainer(SEXP Rv, SEXP RNamed, bool bPrimeFacs, bool bAllFacs, SEXP RNumThreads, int maxThreads);
RcppExport SEXP _RcppAlgos_PollardRhoContainer(SEXP RvSEXP, SEXP RNamedSEXP, SEXP bPrimeFacsSEXP, SEXP bAllFacsSEXP, SEXP RNumThreadsSEXP, SEXP maxThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNamed(RNamedSEXP);
    Rcpp::traits::input_parameter< bool >::type bPrimeFacs(bPrimeFacsSEXP);
    Rcpp::traits::input_parameter< bool >::type bAllFacs(bAllFacsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumThreads(RNumThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxThreads(maxThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(PollardRhoContainer(Rv, RNamed, bPrimeFacs, bAllFacs, RNumThreads, maxThreads));
    return rcpp_result_gen;
END_RCPP
}
// PrimeCountRcpp
SEXP PrimeCountRcpp(SEXP Rn, SEXP RNumThreads, int maxThreads);
RcppExport SEXP _RcppAlgos_PrimeCountRcpp(SEXP RnSEXP, SEXP RNumThreadsSEXP, SEXP maxThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rn(RnSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumThreads(RNumThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxThreads(maxThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(PrimeCountRcpp(Rn, RNumThreads, maxThreads));
    return rcpp_result_gen;
END_RCPP
}
// EratosthenesRcpp
SEXP EratosthenesRcpp(SEXP Rb1, SEXP Rb2, SEXP RNumThreads, int maxCores, int maxThreads);
RcppExport SEXP _RcppAlgos_EratosthenesRcpp(SEXP Rb1SEXP, SEXP Rb2SEXP, SEXP RNumThreadsSEXP, SEXP maxCoresSEXP, SEXP maxThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rb1(Rb1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rb2(Rb2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumThreads(RNumThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxCores(maxCoresSEXP);
    Rcpp::traits::input_parameter< int >::type maxThreads(maxThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(EratosthenesRcpp(Rb1, Rb2, RNumThreads, maxCores, maxThreads));
    return rcpp_result_gen;
END_RCPP
}
// SampleRcpp
SEXP SampleRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, SEXP RFreqs, SEXP RindexVec, bool IsComb, bool IsFactor, SEXP RmySeed, SEXP RNumSamp, Rcpp::Function baseSample, SEXP stdFun, SEXP myEnv, SEXP Rparallel, SEXP RNumThreads, int maxThreads, SEXP RNamed);
RcppExport SEXP _RcppAlgos_SampleRcpp(SEXP RvSEXP, SEXP RmSEXP, SEXP RrepetitionSEXP, SEXP RFreqsSEXP, SEXP RindexVecSEXP, SEXP IsCombSEXP, SEXP IsFactorSEXP, SEXP RmySeedSEXP, SEXP RNumSampSEXP, SEXP baseSampleSEXP, SEXP stdFunSEXP, SEXP myEnvSEXP, SEXP RparallelSEXP, SEXP RNumThreadsSEXP, SEXP maxThreadsSEXP, SEXP RNamedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type Rv(RvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rm(RmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rrepetition(RrepetitionSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RFreqs(RFreqsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RindexVec(RindexVecSEXP);
    Rcpp::traits::input_parameter< bool >::type IsComb(IsCombSEXP);
    Rcpp::traits::input_parameter< bool >::type IsFactor(IsFactorSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RmySeed(RmySeedSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumSamp(RNumSampSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type baseSample(baseSampleSEXP);
    Rcpp::traits::input_parameter< SEXP >::type stdFun(stdFunSEXP);
    Rcpp::traits::input_parameter< SEXP >::type myEnv(myEnvSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Rparallel(RparallelSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNumThreads(RNumThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type maxThreads(maxThreadsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RNamed(RNamedSEXP);
    rcpp_result_gen = Rcpp::wrap(SampleRcpp(Rv, Rm, Rrepetition, RFreqs, RindexVec, IsComb, IsFactor, RmySeed, RNumSamp, baseSample, stdFun, myEnv, Rparallel, RNumThreads, maxThreads, RNamed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RcppAlgos_CheckReturn", (DL_FUNC) &_RcppAlgos_CheckReturn, 7},
    {"_RcppAlgos_CombinatoricsStndrd", (DL_FUNC) &_RcppAlgos_CombinatoricsStndrd, 11},
    {"_RcppAlgos_CombinatoricsApply", (DL_FUNC) &_RcppAlgos_CombinatoricsApply, 9},
    {"_RcppAlgos_CombinatoricsCount", (DL_FUNC) &_RcppAlgos_CombinatoricsCount, 5},
    {"_RcppAlgos_ComboGroupsCountCpp", (DL_FUNC) &_RcppAlgos_ComboGroupsCountCpp, 2},
    {"_RcppAlgos_ComboGroupsRcpp", (DL_FUNC) &_RcppAlgos_ComboGroupsRcpp, 15},
    {"_RcppAlgos_CombinatoricsCnstrt", (DL_FUNC) &_RcppAlgos_CombinatoricsCnstrt, 15},
    {"_RcppAlgos_DivNumSieve", (DL_FUNC) &_RcppAlgos_DivNumSieve, 6},
    {"_RcppAlgos_GetClassVals", (DL_FUNC) &_RcppAlgos_GetClassVals, 8},
    {"_RcppAlgos_MotleyContainer", (DL_FUNC) &_RcppAlgos_MotleyContainer, 6},
    {"_RcppAlgos_cpp11GetNumThreads", (DL_FUNC) &_RcppAlgos_cpp11GetNumThreads, 0},
    {"_RcppAlgos_PollardRhoContainer", (DL_FUNC) &_RcppAlgos_PollardRhoContainer, 6},
    {"_RcppAlgos_PrimeCountRcpp", (DL_FUNC) &_RcppAlgos_PrimeCountRcpp, 3},
    {"_RcppAlgos_EratosthenesRcpp", (DL_FUNC) &_RcppAlgos_EratosthenesRcpp, 5},
    {"_RcppAlgos_SampleRcpp", (DL_FUNC) &_RcppAlgos_SampleRcpp, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_RcppAlgos(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
