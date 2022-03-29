#include "Permutations/PermuteDistinctApply.h"
#include "Permutations/PermuteMultisetApply.h"
#include "Combinations/ComboDistinctApply.h"
#include "Combinations/ComboMultisetApply.h"
#include "Permutations/PermuteRepApply.h"
#include "Combinations/ComboRepApply.h"
#include "CleanConvert.h"

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
        if (!Rf_isVector(RFunVal)) cpp11::stop("'FUN.VALUE' must be a vector");
        const int commonLen = Rf_length(RFunVal);

        switch (TYPEOF(RFunVal)) {
            case STRSXP : {
                SEXP res = PROTECT(Rf_allocVector(STRSXP, nRows * commonLen));
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, STRSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case CPLXSXP : {
                SEXP res = PROTECT(Rf_allocVector(CPLXSXP, nRows * commonLen));
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, CPLXSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case RAWSXP : {
                SEXP res = PROTECT(Rf_allocVector(RAWSXP, nRows * commonLen));
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, RAWSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case LGLSXP : {
                SEXP res = PROTECT(Rf_allocVector(LGLSXP, nRows * commonLen));
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, LGLSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case INTSXP : {
                SEXP res = PROTECT(Rf_allocVector(INTSXP, nRows * commonLen));
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, INTSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } case REALSXP : {
                SEXP res = PROTECT(Rf_allocVector(REALSXP, nRows * commonLen));
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, REALSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                UNPROTECT(1);
                return res;
            } default : {
                SEXP res = PROTECT(Rf_allocVector(VECSXP, nRows * commonLen));
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
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
            if (IsMult) {
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
        if (!Rf_isVector(RFunVal)) cpp11::stop("'FUN.VALUE' must be a vector");
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
            if (IsMult) {
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

SEXP GetCombPermApply(SEXP Rv, const std::vector<double> &vNum,
                      const std::vector<int> &vInt, int n, int m,
                      bool IsComb, bool IsRep, bool IsMult,
                      const std::vector<int> &freqs, std::vector<int> &z,
                      const std::vector<int> &myReps, VecType myType,
                      int nRows, SEXP stdFun, SEXP myEnv, SEXP RFunVal) {

    switch (myType) {
        case VecType::Character : {
            SEXP charVec = PROTECT(Rf_duplicate(Rv));
            SEXP vectorPass = PROTECT(Rf_allocVector(STRSXP, m));

            SEXP res = PROTECT(ApplyFunction(charVec, vectorPass, n, m,
                                             IsComb, IsRep, nRows, freqs,
                                             z, IsMult, stdFun, myEnv,
                                             RFunVal));
            UNPROTECT(3);
            return res;
        } case VecType::Complex : {
            SEXP vectorPass = PROTECT(Rf_allocVector(CPLXSXP, m));
            Rcomplex* ptr_vec = COMPLEX(vectorPass);

            Rcomplex* cmplxVec = COMPLEX(Rv);
            std::vector<Rcomplex> vCmplx(cmplxVec, cmplxVec + n);
            SEXP res = PROTECT(ApplyFunction(vCmplx, vectorPass, ptr_vec,
                                             n, m, IsComb, IsRep, nRows,
                                             freqs, z, IsMult, stdFun,
                                             myEnv, RFunVal));
            UNPROTECT(2);
            return res;
        } case VecType::Raw : {
            SEXP vectorPass = PROTECT(Rf_allocVector(RAWSXP, m));
            Rbyte* ptr_vec = RAW(vectorPass);

            Rbyte* rawVec = RAW(Rv);
            std::vector<Rbyte> vByte(rawVec, rawVec + n);
            SEXP res = PROTECT(ApplyFunction(vByte, vectorPass, ptr_vec,
                                             n, m, IsComb, IsRep, nRows,
                                             freqs, z, IsMult, stdFun,
                                             myEnv, RFunVal));
            UNPROTECT(2);
            return res;
        } case VecType::Logical : {
            SEXP vectorPass = PROTECT(Rf_allocVector(LGLSXP, m));
            int* ptr_vec = LOGICAL(vectorPass);
            SEXP res = PROTECT(ApplyFunction(vInt, vectorPass, ptr_vec,
                                             n, m, IsComb, IsRep, nRows,
                                             freqs, z, IsMult, stdFun,
                                             myEnv, RFunVal));
            UNPROTECT(2);
            return res;
        } case VecType::Integer : {
            SEXP vectorPass = PROTECT(Rf_allocVector(INTSXP, m));
            int* ptr_vec = INTEGER(vectorPass);
            SEXP res = PROTECT(ApplyFunction(vInt, vectorPass, ptr_vec,
                                             n, m, IsComb, IsRep, nRows,
                                             freqs, z, IsMult, stdFun,
                                             myEnv, RFunVal));
            UNPROTECT(2);
            return res;
        } default : {
            SEXP vectorPass = PROTECT(Rf_allocVector(REALSXP, m));
            double* ptr_vec = REAL(vectorPass);
            SEXP res = PROTECT(ApplyFunction(vNum, vectorPass, ptr_vec,
                                             n, m, IsComb, IsRep, nRows,
                                             freqs, z, IsMult, stdFun,
                                             myEnv, RFunVal));
            UNPROTECT(2);
            return res;
        }
    }
}
