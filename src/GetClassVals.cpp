#include "ImportExportMPZ.h"
#include "GetClassValues.h"
#include "ComputedCount.h"
#include "CheckReturn.h"
#include "SetUpUtils.h"

SEXP CopyRv(SEXP Rv, const std::vector<int> &vInt,
            const std::vector<double> &vNum,
            VecType myType, bool IsFactor) {

    if (myType > VecType::Numeric || IsFactor) {
        return Rf_duplicate(Rv);
    } else if (myType == VecType::Integer) {
        return GetIntVec(vInt);
    } else {
        return GetDblVec(vNum);
    }
}

SEXP GetClassVals(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                  SEXP RIsComb, SEXP stdFun, SEXP RThreads,
                  SEXP RmaxThreads, SEXP RIsCnstrd) {

    int n = 0;
    int m = 0;

    bool IsMult = false;
    VecType myType = VecType::Integer;

    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;
    std::vector<double> vNum;

    bool IsRep = CleanConvert::convertFlag(RisRep, "repetition");
    const bool IsComb = CleanConvert::convertFlag(RIsComb, "IsComb");
    const bool IsConstrained = Rf_asLogical(RIsCnstrd);
    const bool IsFactor = Rf_isFactor(Rv);

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum, Rv,
              RFreqs, Rm, n, m, IsMult, IsRep, IsConstrained);

    const SEXP sexpVec = PROTECT(CopyRv(Rv, vInt, vNum, myType, IsFactor));
    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp = (computedRows > SampleLimit);

    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }

    SEXP sexpNumRows = PROTECT(CleanConvert::GetCount(
        IsGmp, computedRowsMpz, computedRows
    ));

    SEXP freqsInfo = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(freqsInfo, 0, GetIntVec(myReps));
    SET_VECTOR_ELT(freqsInfo, 1, GetIntVec(freqs));

    // Needed to determine if nextFullPerm or nextPerm will be called
    const bool IsFullPerm = (IsComb || IsRep) ? false :
        (m == n || m == static_cast<int>(freqs.size()));

    SEXP bVec = PROTECT(Rf_allocVector(LGLSXP, 6));
    INTEGER(bVec)[0] = IsFactor;
    INTEGER(bVec)[1] = IsComb;
    INTEGER(bVec)[2] = IsMult;
    INTEGER(bVec)[3] = IsRep;
    INTEGER(bVec)[4] = IsGmp;
    INTEGER(bVec)[5] = IsFullPerm;

    const bool applyFun = !Rf_isNull(stdFun) && !IsFactor;

    if (applyFun && !Rf_isFunction(stdFun)) {
        Rf_error("FUN must be a function!");
    }

    // RVals is a list containing: v, vNum, vInt, m,
    // RcompRows, nThreads, maxThreads
    SEXP RVals = PROTECT(Rf_allocVector(VECSXP, 7));
    SET_VECTOR_ELT(RVals, 0, sexpVec);
    SET_VECTOR_ELT(RVals, 1, GetDblVec(vNum));
    SET_VECTOR_ELT(RVals, 2, GetIntVec(vInt));
    SET_VECTOR_ELT(RVals, 3, Rf_ScalarInteger(m));
    SET_VECTOR_ELT(RVals, 4, sexpNumRows);
    SET_VECTOR_ELT(RVals, 5, RmaxThreads);
    SET_VECTOR_ELT(RVals, 6, RThreads);

    const char *names[] = {"RVals", "bVec", "FreqsInfo",
                           "applyFun", ""};
    SEXP res = PROTECT(Rf_mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, RVals);
    SET_VECTOR_ELT(res, 1, bVec);
    SET_VECTOR_ELT(res, 2, freqsInfo);
    SET_VECTOR_ELT(res, 3, Rf_ScalarLogical(applyFun));

    UNPROTECT(6);
    return res;
}
