#include "Permutations/PermuteDistinctApply.h"
#include "Permutations/PermuteMultisetApply.h"
#include "Combinations/ComboDistinctApply.h"
#include "Combinations/ComboMultisetApply.h"
#include "Permutations/PermuteRepApply.h"
#include "Combinations/ComboRepApply.h"
#include "CombinatoricsApply.h"
#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "SetUpUtils.h"

template <typename T>
void VecApply(SEXP res, const std::vector<T> &v, SEXP vectorPass,
              T* ptr_vec, int n, int m, bool IsComb, bool IsRep, int nRows,
              const std::vector<int> &freqs, std::vector<int> &z, bool IsMult,
              SEXP stdFun, SEXP rho, int commonLen, int commonType) {

    SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));

    if (IsComb) {
        if (IsMult) {
            MultisetComboApplyFun(res, v, vectorPass, ptr_vec, z, n, m, nRows,
                                  sexpFun, rho, freqs, commonLen, commonType);
        } else if (IsRep) {
            ComboRepApplyFun(res, v, vectorPass, ptr_vec, z, n, m,
                             nRows, sexpFun, rho, commonLen, commonType);
        } else {
            ComboDistinctApplyFun(res, v, vectorPass, ptr_vec, z, n, m,
                                  nRows, sexpFun, rho, commonLen, commonType);
        }
    } else {
        if (IsMult) {
            MultisetPermuteApplyFun(res, v, vectorPass, ptr_vec, z,
                                    n, m, nRows, sexpFun, rho, freqs,
                                    commonLen, commonType);
        } else if (IsRep) {
            PermuteRepApplyFun(res, v, vectorPass, ptr_vec, z, n, m,
                               nRows, sexpFun, rho, commonLen, commonType);
        } else {
            PermuteDistinctApplyFun(res, v, vectorPass, ptr_vec,
                                    z, n, m, nRows, sexpFun, rho,
                                    commonLen, commonType);
        }
    }

    UNPROTECT(1);
}

void VecApply(SEXP res, SEXP v, SEXP vectorPass,
              int n, int m, bool IsComb, bool IsRep, int nRows,
              const std::vector<int> &freqs, std::vector<int> &z, bool IsMult,
              SEXP stdFun, SEXP rho, int commonLen, int commonType) {

    SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));

    if (IsComb) {
        if (IsMult) {
            MultisetComboApplyFun(res, v, vectorPass, z, n, m, nRows,
                                  sexpFun, rho, freqs, commonLen, commonType);
        } else if (IsRep) {
            ComboRepApplyFun(res, v, vectorPass, z, n, m, nRows,
                             sexpFun, rho, commonLen, commonType);
        } else {
            ComboDistinctApplyFun(res, v, vectorPass, z, n, m, nRows,
                                  sexpFun, rho, commonLen, commonType);
        }
    } else {
        if (IsMult) {
            MultisetPermuteApplyFun(res, v, vectorPass, z, n, m,
                                    nRows, sexpFun, rho, freqs,
                                    commonLen, commonType);
        } else if (IsRep) {
            PermuteRepApplyFun(res, v, vectorPass, z, n, m, nRows,
                               sexpFun, rho, commonLen, commonType);
        } else {
            PermuteDistinctApplyFun(res, v, vectorPass, z,
                                    n, m, nRows, sexpFun, rho,
                                    commonLen, commonType);
        }
    }

    UNPROTECT(1);
}

SEXP ApplyFunction(SEXP v, SEXP vectorPass, int n, int m, bool IsComb,
                   bool IsRep, int nRows, const std::vector<int> &freqs,
                   std::vector<int> &z, bool IsMult, SEXP stdFun,
                   SEXP rho, SEXP RFunVal) {

    if (!Rf_isNull(RFunVal)) {
        if (!Rf_isVector(RFunVal)) Rf_error("'FUN.VALUE' must be a vector");
        const int commonLen = Rf_length(RFunVal);

        switch (TYPEOF(RFunVal)) {
            case STRSXP : {
                SEXP res = PROTECT(Rf_allocVector(STRSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, STRSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case CPLXSXP : {
                SEXP res = PROTECT(Rf_allocVector(CPLXSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, CPLXSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case RAWSXP : {
                SEXP res = PROTECT(Rf_allocVector(RAWSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, RAWSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case LGLSXP : {
                SEXP res = PROTECT(Rf_allocVector(LGLSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, LGLSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case INTSXP : {
                SEXP res = PROTECT(Rf_allocVector(INTSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, INTSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case REALSXP : {
                SEXP res = PROTECT(Rf_allocVector(REALSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, REALSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } default : {
                SEXP res = PROTECT(Rf_allocVector(VECSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, VECSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            }
        }
    } else {
        SEXP myList = PROTECT(Rf_allocVector(VECSXP, nRows));
        SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));

        if (IsComb) {
            if (IsMult){
                MultisetComboApplyFun(myList, v, vectorPass, z, n,
                                      m, nRows, sexpFun, rho, freqs);
            } else if (IsRep) {
                ComboRepApplyFun(myList, v, vectorPass, z,
                                 n, m, nRows, sexpFun, rho);
            } else {
                ComboDistinctApplyFun(myList, v, vectorPass, z,
                                      n, m, nRows, sexpFun, rho);
            }
        } else {
            if (IsMult) {
                MultisetPermuteApplyFun(myList, v, vectorPass, z, n,
                                        m, nRows, sexpFun, rho, freqs);
            } else if (IsRep) {
                PermuteRepApplyFun(myList, v, vectorPass, z,
                                   n, m, nRows, sexpFun, rho);
            } else {
                PermuteDistinctApplyFun(myList, v, vectorPass, z,
                                        n, m, nRows, sexpFun, rho);
            }
        }

        UNPROTECT(2);
        return myList;
    }
}

template <typename T>
SEXP ApplyFunction(const std::vector<T> &v, SEXP vectorPass,
                   T* ptr_vec, int n, int m, bool IsComb, bool IsRep,
                   int nRows, const std::vector<int> &freqs,
                   std::vector<int> &z, bool IsMult, SEXP stdFun,
                   SEXP rho, SEXP RFunVal) {

    if (!Rf_isNull(RFunVal)) {
        if (!Rf_isVector(RFunVal)) Rf_error("'FUN.VALUE' must be a vector");
        const int commonLen = Rf_length(RFunVal);

        switch (TYPEOF(RFunVal)) {
            case STRSXP : {
                SEXP res = PROTECT(Rf_allocVector(STRSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, STRSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case CPLXSXP : {
                SEXP res = PROTECT(Rf_allocVector(CPLXSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, CPLXSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case RAWSXP : {
                SEXP res = PROTECT(Rf_allocVector(RAWSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, RAWSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case LGLSXP : {
                SEXP res = PROTECT(Rf_allocVector(LGLSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, LGLSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case INTSXP : {
                SEXP res = PROTECT(Rf_allocVector(INTSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, INTSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case REALSXP : {
                SEXP res = PROTECT(Rf_allocVector(REALSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, REALSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } default : {
                SEXP res = PROTECT(Rf_allocVector(VECSXP, nRows * commonLen));

                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, VECSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            }
        }
    } else {
        SEXP myList = PROTECT(Rf_allocVector(VECSXP, nRows));
        SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));

        if (IsComb) {
            if (IsMult){
                MultisetComboApplyFun(myList, v, vectorPass, ptr_vec, z,
                                      n, m, nRows, sexpFun, rho, freqs);
            } else if (IsRep) {
                ComboRepApplyFun(myList, v, vectorPass, ptr_vec,
                                 z, n, m, nRows, sexpFun, rho);
            } else {
                ComboDistinctApplyFun(myList, v, vectorPass, ptr_vec,
                                      z, n, m, nRows, sexpFun, rho);
            }
        } else {
            if (IsMult) {
                MultisetPermuteApplyFun(myList, v, vectorPass, ptr_vec, z,
                                        n, m, nRows, sexpFun, rho, freqs);
            } else if (IsRep) {
                PermuteRepApplyFun(myList, v, vectorPass, ptr_vec,
                                   z, n, m, nRows, sexpFun, rho);
            } else {
                PermuteDistinctApplyFun(myList, v, vectorPass, ptr_vec,
                                        z, n, m, nRows, sexpFun, rho);
            }
        }

        UNPROTECT(2);
        return myList;
    }
}

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

    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    const bool IsComb = CleanConvert::convertLogical(RIsComb, "IsComb");

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
    SetNumResults(IsGmp, bLower, bUpper, false, upperMpz.get(), lowerMpz.get(),
                  lower, upper, computedRows, computedRowsMpz, nRows, userNumRows);

    switch (myType) {
        case VecType::Character : {
            SEXP vectorPass = PROTECT(Rf_allocVector(STRSXP, m));
            SEXP res = ApplyFunction(Rv, vectorPass, n, m, IsComb, IsRep,
                                     nRows, freqs, startZ, IsMult, stdFun,
                                     myEnv, RFunVal);
            UNPROTECT(1);
            return res;
        } case VecType::Complex : {
            SEXP vectorPass = PROTECT(Rf_allocVector(CPLXSXP, m));
            Rcomplex* ptr_vec = COMPLEX(vectorPass);

            Rcomplex* cmplxVec = COMPLEX(Rv);
            std::vector<Rcomplex> vCmplx(cmplxVec, cmplxVec + n);
            SEXP res = ApplyFunction(vCmplx, vectorPass, ptr_vec, n, m,
                                     IsComb, IsRep, nRows, freqs, startZ,
                                     IsMult, stdFun, myEnv, RFunVal);
            UNPROTECT(1);
            return res;
        } case VecType::Raw : {
            SEXP vectorPass = PROTECT(Rf_allocVector(RAWSXP, m));
            Rbyte* ptr_vec = RAW(vectorPass);

            Rbyte* rawVec = RAW(Rv);
            std::vector<Rbyte> vByte(rawVec, rawVec + n);
            SEXP res = ApplyFunction(vByte, vectorPass, ptr_vec, n, m,
                                     IsComb, IsRep, nRows, freqs, startZ,
                                     IsMult, stdFun, myEnv, RFunVal);
            UNPROTECT(1);
            return res;
        } case VecType::Logical : {
            SEXP vectorPass = PROTECT(Rf_allocVector(LGLSXP, m));
            int* ptr_vec = LOGICAL(vectorPass);
            SEXP res = ApplyFunction(vInt, vectorPass, ptr_vec, n, m,
                                     IsComb, IsRep, nRows, freqs, startZ,
                                     IsMult, stdFun, myEnv, RFunVal);
            UNPROTECT(1);
            return res;
        } case VecType::Integer : {
            SEXP vectorPass = PROTECT(Rf_allocVector(INTSXP, m));
            int* ptr_vec = INTEGER(vectorPass);
            SEXP res = ApplyFunction(vInt, vectorPass, ptr_vec, n, m,
                                     IsComb, IsRep, nRows, freqs, startZ,
                                     IsMult, stdFun, myEnv, RFunVal);
            UNPROTECT(1);
            return res;
        } default : {
            SEXP vectorPass = PROTECT(Rf_allocVector(REALSXP, m));
            double* ptr_vec = REAL(vectorPass);
            SEXP res = ApplyFunction(vNum, vectorPass, ptr_vec, n, m,
                                     IsComb, IsRep, nRows, freqs, startZ,
                                     IsMult, stdFun, myEnv, RFunVal);
            UNPROTECT(1);
            return res;
        }
    }
}
