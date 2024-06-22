#include "Permutations/PermuteDistinctApply.h"
#include "Permutations/PermuteMultisetApply.h"
#include "Combinations/ComboDistinctApply.h"
#include "Combinations/ComboMultisetApply.h"
#include "Permutations/PermuteRepApply.h"
#include "Combinations/ComboRepApply.h"
#include "CppConvert.h"

template <typename T>
void VecApply(SEXP res, const std::vector<T> &v, SEXP vectorPass,
              T* ptr_vec, int n, int m, bool IsComb, bool IsRep, int nRows,
              const std::vector<int> &freqs, std::vector<int> &z, bool IsMult,
              SEXP stdFun, SEXP rho, int commonLen, int commonType) {

    cpp11::sexp sexpFun = Rf_lang2(stdFun, R_NilValue);
    MARK_NOT_MUTABLE(sexpFun);

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

}

void VecApply(SEXP res, SEXP v, SEXP vectorPass,
              int n, int m, bool IsComb, bool IsRep, int nRows,
              const std::vector<int> &freqs, std::vector<int> &z, bool IsMult,
              SEXP stdFun, SEXP rho, int commonLen, int commonType) {

    cpp11::sexp sexpFun = Rf_lang2(stdFun, R_NilValue);
    MARK_NOT_MUTABLE(sexpFun);

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
                cpp11::sexp res = Rf_allocVector(STRSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, STRSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case CPLXSXP : {
                cpp11::sexp res = Rf_allocVector(CPLXSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, CPLXSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case RAWSXP : {
                cpp11::sexp res = Rf_allocVector(RAWSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, RAWSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case LGLSXP : {
                cpp11::sexp res = Rf_allocVector(LGLSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, LGLSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case INTSXP : {
                cpp11::sexp res = Rf_allocVector(INTSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, INTSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case REALSXP : {
                cpp11::sexp res = Rf_allocVector(REALSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, REALSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } default : {
                cpp11::sexp res = Rf_allocVector(VECSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, n, m, IsComb,
                         IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, VECSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            }
        }
    } else {
        cpp11::sexp myList = Rf_allocVector(VECSXP, nRows);
        cpp11::sexp sexpFun = Rf_lang2(stdFun, R_NilValue);
        MARK_NOT_MUTABLE(sexpFun);

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
                cpp11::sexp res = Rf_allocVector(STRSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, STRSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case CPLXSXP : {
                cpp11::sexp res = Rf_allocVector(CPLXSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, CPLXSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case RAWSXP : {
                cpp11::sexp res = Rf_allocVector(RAWSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, RAWSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case LGLSXP : {
                cpp11::sexp res = Rf_allocVector(LGLSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, LGLSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case INTSXP : {
                cpp11::sexp res = Rf_allocVector(INTSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, INTSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case REALSXP : {
                cpp11::sexp res = Rf_allocVector(REALSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, REALSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } default : {
                cpp11::sexp res = Rf_allocVector(VECSXP, nRows * commonLen);
                VecApply(res, v, vectorPass, ptr_vec, n, m,
                         IsComb, IsRep, nRows, freqs, z, IsMult,
                         stdFun, rho, commonLen, VECSXP);
                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            }
        }
    } else {
        cpp11::sexp myList = Rf_allocVector(VECSXP, nRows);
        cpp11::sexp sexpFun = Rf_lang2(stdFun, R_NilValue);
        MARK_NOT_MUTABLE(sexpFun);

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
            cpp11::sexp charVec = Rf_duplicate(Rv);
            cpp11::sexp vectorPass = Rf_allocVector(STRSXP, m);
            cpp11::sexp res = ApplyFunction(charVec, vectorPass, n, m,
                                            IsComb, IsRep, nRows, freqs,
                                            z, IsMult, stdFun, myEnv,
                                            RFunVal);
            return res;
        } case VecType::Complex : {
            cpp11::sexp vectorPass = Rf_allocVector(CPLXSXP, m);
            Rcomplex* ptr_vec = COMPLEX(vectorPass);
            std::vector<Rcomplex> vCmplx = CppConvert::GetVec<Rcomplex>(Rv);
            cpp11::sexp res = ApplyFunction(vCmplx, vectorPass, ptr_vec,
                                            n, m, IsComb, IsRep, nRows,
                                            freqs, z, IsMult, stdFun,
                                            myEnv, RFunVal);
            return res;
        } case VecType::Raw : {
            cpp11::sexp vectorPass = Rf_allocVector(RAWSXP, m);
            Rbyte* ptr_vec = RAW(vectorPass);
            std::vector<Rbyte> vByte = CppConvert::GetVec<Rbyte>(Rv);
            cpp11::sexp res = ApplyFunction(vByte, vectorPass, ptr_vec,
                                            n, m, IsComb, IsRep, nRows,
                                            freqs, z, IsMult, stdFun,
                                            myEnv, RFunVal);
            return res;
        } case VecType::Logical : {
            cpp11::sexp vectorPass = Rf_allocVector(LGLSXP, m);
            int* ptr_vec = LOGICAL(vectorPass);
            cpp11::sexp res = ApplyFunction(vInt, vectorPass, ptr_vec,
                                            n, m, IsComb, IsRep, nRows,
                                            freqs, z, IsMult, stdFun,
                                            myEnv, RFunVal);
            return res;
        } case VecType::Integer : {
            cpp11::sexp vectorPass = Rf_allocVector(INTSXP, m);
            int* ptr_vec = INTEGER(vectorPass);
            cpp11::sexp res = ApplyFunction(vInt, vectorPass, ptr_vec,
                                            n, m, IsComb, IsRep, nRows,
                                            freqs, z, IsMult, stdFun,
                                            myEnv, RFunVal);
            return res;
        } default : {
            cpp11::sexp vectorPass = Rf_allocVector(REALSXP, m);
            double* ptr_vec = REAL(vectorPass);
            cpp11::sexp res = ApplyFunction(vNum, vectorPass, ptr_vec,
                                            n, m, IsComb, IsRep, nRows,
                                            freqs, z, IsMult, stdFun,
                                            myEnv, RFunVal);
            return res;
        }
    }
}
