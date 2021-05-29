#include "CombinatoricsGlue.h"
#include "CombinatoricsMain.h"
#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "SetUpUtils.h"

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

    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    const bool IsComb = CleanConvert::convertLogical(RIsComb, "IsComb");
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
    SetNumResults(IsGmp, bLower, bUpper, false, upperMpz.get(),
                  lowerMpz.get(), lower, upper, computedRows,
                  computedRowsMpz, nRows, userNumRows);

    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows,
               myType, nThreads, RNumThreads, limit);

    int phaseOne = 0;
    bool generalRet = true;
    const bool IsCharacter = myType == VecType::Character;

    PermuteSpecific(phaseOne, generalRet, n, m, nRows,
                    IsMult, IsCharacter, IsComb, bLower, IsRep);

    switch (myType) {
        case VecType::Character : {
            SEXP res = PROTECT(Rf_allocMatrix(STRSXP, nRows, m));

            CharacterGlue(res, Rv, IsComb, startZ, n, m,
                          nRows, freqs, IsMult, IsRep);

            UNPROTECT(1);
            return res;
        } case VecType::Complex : {
            std::vector<Rcomplex> stlCmplxVec(n);
            Rcomplex* vecCmplx = COMPLEX(Rv);

            for (int i = 0; i < n; ++i)
                stlCmplxVec[i] = vecCmplx[i];

            SEXP res = PROTECT(Rf_allocMatrix(CPLXSXP, nRows, m));
            Rcomplex* matCmplx = COMPLEX(res);

            ManagerGlue(matCmplx, stlCmplxVec, startZ, n, m, nRows,
                        IsComb, phaseOne, generalRet, freqs, IsMult, IsRep);

            UNPROTECT(1);
            return res;
        }
        case VecType::Raw : {
            std::vector<Rbyte> stlRawVec(n);
            Rbyte* rawVec = RAW(Rv);

            for (int i = 0; i < n; ++i)
                stlRawVec[i] = rawVec[i];

            SEXP res = PROTECT(Rf_allocMatrix(RAWSXP, nRows, m));
            Rbyte* rawMat = RAW(res);

            ManagerGlue(rawMat, stlRawVec, startZ, n, m, nRows,
                        IsComb, phaseOne, generalRet, freqs, IsMult, IsRep);

            UNPROTECT(1);
            return res;
        } case VecType::Logical : {
            vInt.assign(n, 0);
            int* vecBool = LOGICAL(Rv);

            for (int i = 0; i < n; ++i)
                vInt[i] = vecBool[i];

            SEXP res = PROTECT(Rf_allocMatrix(LGLSXP, nRows, m));
            int* matBool = LOGICAL(res);

            ManagerGlue(matBool, vInt, startZ, n, m, nRows, IsComb,
                        phaseOne, generalRet, freqs, IsMult, IsRep);

            UNPROTECT(1);
            return res;
        } case VecType::Integer : {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, nRows, m));
            int* matInt = INTEGER(res);

            ParallelGlue(matInt, vInt, n, m, phaseOne, generalRet, IsComb,
                         Parallel, IsRep, IsMult, IsGmp, freqs, startZ,
                         myReps, lower, lowerMpz[0], nRows, nThreads);

            if (Rf_isFactor(Rv)) {
                SetFactorClass(res, Rv);
            }

            UNPROTECT(1);
            return res;
        } default : {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, nRows, m));
            double* matNum = REAL(res);

            ParallelGlue(matNum, vNum, n, m, phaseOne, generalRet, IsComb,
                         Parallel, IsRep, IsMult, IsGmp, freqs, startZ,
                         myReps, lower, lowerMpz[0], nRows, nThreads);

            UNPROTECT(1);
            return res;
        }
    }
}
