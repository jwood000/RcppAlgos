#include "GetCombPermApply.h"
#include "ComputedCount.h"
#include "SetUpUtils.h"

[[cpp11::register]]
SEXP CombinatoricsApply(SEXP Rv, SEXP Rm, SEXP RisRep,
                        SEXP RFreqs, SEXP Rlow, SEXP Rhigh,
                        SEXP stdFun, SEXP myEnv,
                        SEXP RFunVal, SEXP RIsComb) {
    int n = 0;
    int m = 0;
    int nRows = 0;

    VecType myType = VecType::Integer;
    bool IsMult = false;

    std::vector<double> vNum;
    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;

    bool IsRep = CppConvert::convertFlag(RisRep, "repetition");
    const bool IsComb = CppConvert::convertFlag(RIsComb, "IsComb");

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum,
              Rv, RFreqs, Rm, n, m, IsMult, IsRep);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp = (computedRows > Significand53);

    mpz_class computedRowsMpz;

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult, IsComb,
                          IsRep, n, m, Rm, freqs, myReps);
    }

    double lower = 0;
    double upper = 0;

    bool bLower = false;
    bool bUpper = false;

    mpz_class lowerMpz;
    mpz_class upperMpz;

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz, upperMpz, computedRowsMpz, computedRows);

    std::vector<int> startZ(m);
    SetStartZ(myReps, freqs, startZ, IsComb, n, m,
              lower, lowerMpz, IsRep, IsMult, IsGmp);

    double userNumRows = 0;   // IsGenCnstrd = false
    SetNumResults(IsGmp, bLower, bUpper, true, upperMpz,
                  lowerMpz, lower, upper, computedRows,
                  computedRowsMpz, nRows, userNumRows);

    return GetCombPermApply(Rv, vNum, vInt, n, m, IsComb, IsRep,
                            IsMult, freqs, startZ, myReps, myType,
                            nRows, stdFun, myEnv, RFunVal);
}
