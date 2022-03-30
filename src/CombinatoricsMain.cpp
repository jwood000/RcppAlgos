#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "GetCombPerm.h"
#include "SetUpUtils.h"

[[cpp11::register]]
SEXP CombinatoricsStndrd(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                         SEXP Rlow, SEXP Rhigh, SEXP Rparallel,
                         SEXP RNumThreads, SEXP RmaxThreads, SEXP RIsComb) {
    int n = 0;
    int m = 0;
    int nRows = 0;
    int nThreads = 1;
    int maxThreads = 1;

    VecType myType = VecType::Integer;
    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

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
    const bool IsGmp = (computedRows > Significand53);

    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult, IsComb,
                          IsRep, n, m, Rm, freqs, myReps);
    }

    double lower = 0;
    double upper = 0;

    bool bLower = false;
    bool bUpper = false;

    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);

    mpz_init(lowerMpz[0]);
    mpz_init(upperMpz[0]);

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz.get(), upperMpz.get(), computedRowsMpz, computedRows);

    std::vector<int> startZ(m);
    SetStartZ(myReps, freqs, startZ, IsComb, n, m,
              lower, lowerMpz[0], IsRep, IsMult, IsGmp);

    double userNumRows = 0;
    SetNumResults(IsGmp, bLower, bUpper, true, upperMpz[0],
                  lowerMpz[0], lower, upper, computedRows,
                  computedRowsMpz, nRows, userNumRows);

    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows,
               myType, nThreads, RNumThreads, limit);

    int phaseOne = 0;
    bool generalRet = true;
    const bool IsCharacter = myType == VecType::Character;

    PermuteSpecific(phaseOne, generalRet, n, m, nRows,
                    IsMult, IsCharacter, IsComb, bLower, IsRep);

    return GetCombPerms(Rv, vNum, vInt, n, m, phaseOne, generalRet,
                        IsComb, Parallel, IsRep, IsMult, IsGmp, freqs,
                        startZ, myReps, lower, lowerMpz[0], nRows,
                        nThreads, myType);
}
