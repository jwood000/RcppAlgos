#include "GetCombPermApply.h"
#include "Cpp14MakeUnique.h"
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

    bool IsRep = CleanConvert::convertFlag(RisRep, "repetition");
    const bool IsComb = CleanConvert::convertFlag(RIsComb, "IsComb");

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum,
              Rv, RFreqs, Rm, n, m, IsMult, IsRep);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp = (computedRows > Significand53);

    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult, IsComb,
                          IsRep, n, m, Rm, freqs, myReps);
    }

    double lower = 0, upper = 0;
    bool bLower = false, bUpper = false;

    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);

    mpz_init(lowerMpz[0]);
    mpz_init(upperMpz[0]);

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz.get(), upperMpz.get(), computedRowsMpz, computedRows);

    std::vector<int> startZ(m);
    SetStartZ(myReps, freqs, startZ, IsComb, n, m,
              lower, lowerMpz[0], IsRep, IsMult, IsGmp);

    double userNumRows = 0;   // IsGenCnstrd = false
    SetNumResults(IsGmp, bLower, bUpper, true, upperMpz[0],
                  lowerMpz[0], lower, upper, computedRows,
                  computedRowsMpz, nRows, userNumRows);

    mpz_clear(computedRowsMpz);
    return GetCombPermApply(Rv, vNum, vInt, n, m, IsComb, IsRep,
                            IsMult, freqs, startZ, myReps, myType,
                            nRows, stdFun, myEnv, RFunVal);
}
