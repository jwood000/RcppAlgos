#include "ClassUtils/PrevCombinatorics.h"
#include "SetUpUtils.h"
#include "FunAssign.h"

void VecApplyPrev(SEXP res, SEXP v, SEXP vectorPass,
                  std::vector<int> &z, prevIterPtr prevIter, int n,
                  int m, int nRows, const std::vector<int> &freqs,
                  bool IsComb, bool IsMult, SEXP sexpFun, SEXP rho,
                  int commonLen = 1, int commonType = INTSXP) {

    const int n1 = IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1);
    const int lastRow = nRows - 1;
    const int retType = TYPEOF(res);

    // We iterator to the penultimate row to avoid iterating z one too many times
    for (int count = 0, m1 = m - 1; count < lastRow; ++count,
         prevIter(freqs, z, n1, m1)) {

        for (int j = 0; j < m; ++j) {
            SET_STRING_ELT(vectorPass, j, STRING_ELT(v, z[j]));
        }

        FunAssign(res, vectorPass, sexpFun, rho,
                  commonType, commonLen, count, nRows, retType);
    }

    // Get the last result
    for (int j = 0; j < m; ++j) {
        SET_STRING_ELT(vectorPass, j, STRING_ELT(v, z[j]));
    }

    FunAssign(res, vectorPass, sexpFun, rho,
              commonType, commonLen, lastRow, nRows, retType);
}

template <typename T>
void VecApplyPrev(SEXP res, const std::vector<T> &v,
                  SEXP vectorPass, T* ptr_vec, std::vector<int> &z,
                  prevIterPtr prevIter, int n, int m, int nRows,
                  const std::vector<int> &freqs, bool IsComb,
                  bool IsMult, SEXP sexpFun, SEXP rho,
                  int commonLen = 1, int commonType = INTSXP) {

    const int n1 = IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1);
    const int lastRow = nRows - 1;
    const int retType = TYPEOF(res);

    // We iterater to the pentultimate row to avoid iterating z one too many times
    for (int count = 0, m1 = m - 1; count < lastRow; ++count,
         prevIter(freqs, z, n1, m1)) {

        for (int j = 0; j < m; ++j) {
            ptr_vec[j] = v[z[j]];
        }

        FunAssign(res, vectorPass, sexpFun, rho,
                  commonType, commonLen, count, nRows, retType);
    }

    // Get the last result
    for (int j = 0; j < m; ++j) {
        ptr_vec[j] = v[z[j]];
    }

    FunAssign(res, vectorPass, sexpFun, rho,
              commonType, commonLen, lastRow, nRows, retType);
}

SEXP ApplyFunPrev(SEXP v, SEXP vectorPass, const std::vector<int> &freqs,
                  std::vector<int> &z, SEXP stdFun, SEXP rho, SEXP RFunVal,
                  prevIterPtr prevIter, int n, int m, int nRows,
                  bool IsComb, bool IsMult) {

    cpp11::sexp sexpFun = Rf_lang2(stdFun, R_NilValue);

    if (!Rf_isNull(RFunVal)) {
        if (!Rf_isVector(RFunVal)) cpp11::stop("'FUN.VALUE' must be a vector");
        const int commonLen = Rf_length(RFunVal);

        switch (TYPEOF(RFunVal)) {
            case STRSXP : {
                cpp11::sexp res = Rf_allocVector(STRSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, z, prevIter, n, m,
                             nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, STRSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case CPLXSXP : {
                cpp11::sexp res = Rf_allocVector(CPLXSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, z, prevIter, n, m,
                             nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, CPLXSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case RAWSXP : {
                cpp11::sexp res = Rf_allocVector(RAWSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, z, prevIter, n, m,
                             nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, RAWSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case LGLSXP : {
                cpp11::sexp res = Rf_allocVector(LGLSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, z, prevIter, n, m,
                             nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, LGLSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case INTSXP : {
                cpp11::sexp res = Rf_allocVector(INTSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, z, prevIter, n, m,
                             nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, INTSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case REALSXP : {
                cpp11::sexp res = Rf_allocVector(REALSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, z, prevIter, n, m,
                             nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, REALSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } default : {
                cpp11::sexp res = Rf_allocVector(VECSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, z, prevIter, n, m,
                             nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, VECSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            }
        }
    } else {
        cpp11::sexp myList = Rf_allocVector(VECSXP, nRows);
        VecApplyPrev(myList, v, vectorPass, z, prevIter, n, m,
                     nRows, freqs, IsComb, IsMult, sexpFun, rho);
        return myList;
    }
}

template <typename T>
SEXP ApplyFunPrev(const std::vector<T> &v, SEXP vectorPass, T* ptr_vec,
                  const std::vector<int> &freqs, std::vector<int> &z,
                  SEXP stdFun, SEXP rho, SEXP RFunVal, prevIterPtr prevIter,
                  int n, int m, int nRows, bool IsComb, bool IsMult) {

    cpp11::sexp sexpFun = Rf_lang2(stdFun, R_NilValue);

    if (!Rf_isNull(RFunVal)) {
        if (!Rf_isVector(RFunVal)) cpp11::stop("'FUN.VALUE' must be a vector");
        const int commonLen = Rf_length(RFunVal);

        switch (TYPEOF(RFunVal)) {
            case STRSXP : {
                cpp11::sexp res = Rf_allocVector(STRSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, ptr_vec, z, prevIter,
                             n, m, nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, STRSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case CPLXSXP : {
                cpp11::sexp res = Rf_allocVector(CPLXSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, ptr_vec, z, prevIter,
                             n, m, nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, CPLXSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case RAWSXP : {
                cpp11::sexp res = Rf_allocVector(RAWSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, ptr_vec, z, prevIter,
                             n, m, nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, RAWSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case LGLSXP : {
                cpp11::sexp res = Rf_allocVector(LGLSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, ptr_vec, z, prevIter,
                             n, m, nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, LGLSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case INTSXP : {
                cpp11::sexp res = Rf_allocVector(INTSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, ptr_vec, z, prevIter,
                             n, m, nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, INTSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } case REALSXP : {
                cpp11::sexp res = Rf_allocVector(REALSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, ptr_vec, z, prevIter,
                             n, m, nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, REALSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            } default : {
                cpp11::sexp res = Rf_allocVector(VECSXP, nRows * commonLen);

                VecApplyPrev(res, v, vectorPass, ptr_vec, z, prevIter,
                             n, m, nRows, freqs, IsComb, IsMult, sexpFun,
                             rho, commonLen, VECSXP);

                SetDims(RFunVal, res, commonLen, nRows);
                return res;
            }
        }
    } else {
        cpp11::sexp myList = Rf_allocVector(VECSXP, nRows);
        VecApplyPrev(myList, v, vectorPass, ptr_vec, z, prevIter, n,
                     m, nRows, freqs, IsComb, IsMult, sexpFun, rho);
        return myList;
    }
}

SEXP GetPrevCombPermApply(SEXP Rv, const std::vector<double> &vNum,
                          const std::vector<int> &vInt,
                          const std::vector<int> &myReps,
                          const std::vector<int> &freqs, std::vector<int> &z,
                          prevIterPtr prevIter, int n, int m, bool IsComb,
                          bool IsMult, int nRows, VecType myType,
                          SEXP stdFun, SEXP myEnv, SEXP RFunVal) {

    switch (myType) {
        case VecType::Character : {
            cpp11::sexp charVec = Rf_duplicate(Rv);
            cpp11::sexp vectorPass = Rf_allocVector(STRSXP, m);
            cpp11::sexp res = ApplyFunPrev(charVec, vectorPass, freqs, z,
                                           stdFun, myEnv, RFunVal, prevIter,
                                           n, m, nRows, IsComb, IsMult);
            return res;
        } case VecType::Complex : {
            cpp11::sexp vectorPass = Rf_allocVector(CPLXSXP, m);
            Rcomplex* ptr_vec = COMPLEX(vectorPass);
            std::vector<Rcomplex> vCmplx = CppConvert::GetVec<Rcomplex>(Rv);
            cpp11::sexp res = ApplyFunPrev(vCmplx, vectorPass, ptr_vec,
                                           freqs, z, stdFun, myEnv,
                                           RFunVal, prevIter, n, m,
                                           nRows, IsComb, IsMult);
            return res;
        } case VecType::Raw : {
            cpp11::sexp vectorPass = Rf_allocVector(RAWSXP, m);
            Rbyte* ptr_vec = RAW(vectorPass);
            std::vector<Rbyte> vByte = CppConvert::GetVec<Rbyte>(Rv);
            cpp11::sexp res = ApplyFunPrev(vByte, vectorPass, ptr_vec,
                                           freqs, z, stdFun, myEnv,
                                           RFunVal, prevIter, n, m,
                                           nRows, IsComb, IsMult);
            return res;
        } case VecType::Logical : {
            cpp11::sexp vectorPass = Rf_allocVector(LGLSXP, m);
            int* ptr_vec = LOGICAL(vectorPass);
            cpp11::sexp res = ApplyFunPrev(vInt, vectorPass, ptr_vec,
                                           freqs, z, stdFun, myEnv,
                                           RFunVal, prevIter, n, m,
                                           nRows, IsComb, IsMult);
            return res;
        } case VecType::Integer : {
            cpp11::sexp vectorPass = Rf_allocVector(INTSXP, m);
            int* ptr_vec = INTEGER(vectorPass);
            cpp11::sexp res = ApplyFunPrev(vInt, vectorPass, ptr_vec,
                                           freqs, z, stdFun, myEnv,
                                           RFunVal, prevIter, n, m,
                                           nRows, IsComb, IsMult);
            return res;
        } default : {
            cpp11::sexp vectorPass = Rf_allocVector(REALSXP, m);
            double* ptr_vec = REAL(vectorPass);
            cpp11::sexp res = ApplyFunPrev(vNum, vectorPass, ptr_vec,
                                           freqs, z, stdFun, myEnv,
                                           RFunVal, prevIter, n, m,
                                           nRows, IsComb, IsMult);
            return res;
        }
    }
}
