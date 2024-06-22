#include "SetUpUtils.h"
#include "FunAssign.h"
#include "NthResult.h"

template <typename T>
void SampleApplyFun(SEXP res, const std::vector<T> &v, SEXP vectorPass,
                    T* ptr_vec, const std::vector<double> &mySample,
                    const std::vector<mpz_class> &myBigSamp,
                    const std::vector<int> &myReps,
                    SEXP func, SEXP rho, nthResultPtr nthResFun, int m,
                    int sampSize, bool IsNamed, bool IsGmp, int lenV,
                    int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);
    cpp11::sexp sexpFun = Rf_lang2(func, R_NilValue);

    if (IsGmp) {
        for (int count = 0; count < sampSize; ++count) {
            const std::vector<int> z = nthResFun(lenV, m, 0.0,
                                                 myBigSamp[count], myReps);

            for (int j = 0; j < m; ++j) {
                ptr_vec[j] = v[z[j]];
            }

            FunAssign(res, vectorPass, sexpFun, rho, commonType,
                      commonLen, count, sampSize, retType);
        }
    } else {
        mpz_class mpzDefault;

        for (int count = 0; count < sampSize; ++count) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[count],
                                                 mpzDefault, myReps);

            for (int j = 0; j < m; ++j) {
                ptr_vec[j] = v[z[j]];
            }

            FunAssign(res, vectorPass, sexpFun, rho, commonType,
                      commonLen, count, sampSize, retType);
        }
    }

    SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp, IsNamed);
}

void SampleApplyFun(SEXP res, SEXP v, SEXP vectorPass,
                    const std::vector<double> &mySample,
                    const std::vector<mpz_class> &myBigSamp,
                    const std::vector<int> &myReps,
                    SEXP func, SEXP rho, nthResultPtr nthResFun, int m,
                    int sampSize, bool IsNamed, bool IsGmp, int lenV,
                    int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);
    cpp11::sexp sexpFun = Rf_lang2(func, R_NilValue);

    if (IsGmp) {
        for (int count = 0; count < sampSize; ++count) {
            const std::vector<int> z = nthResFun(lenV, m, 0.0,
                                                 myBigSamp[count], myReps);

            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(vectorPass, j, STRING_ELT(v, z[j]));
            }

            FunAssign(res, vectorPass, sexpFun, rho, commonType,
                      commonLen, count, sampSize, retType);
        }
    } else {
        mpz_class mpzDefault;

        for (int count = 0; count < sampSize; ++count) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[count],
                                                 mpzDefault, myReps);

            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(vectorPass, j, STRING_ELT(v, z[j]));
            }

            FunAssign(res, vectorPass, sexpFun, rho, commonType,
                      commonLen, count, sampSize, retType);
        }
    }

    SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp, IsNamed);
}

template <typename T>
SEXP ApplyFunction(const std::vector<T> &v, SEXP vectorPass,
                   T* ptr_vec, const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   const std::vector<int> &myReps,
                   SEXP stdFun, SEXP rho, SEXP RFunVal,
                   nthResultPtr nthResFun, int m, int n_samp,
                   bool IsNamed, bool IsGmp, int lenV) {

    if (!Rf_isNull(RFunVal)) {
        if (!Rf_isVector(RFunVal)) cpp11::stop("'FUN.VALUE' must be a vector");
        const int commonLen = Rf_length(RFunVal);

        switch (TYPEOF(RFunVal)) {
            case STRSXP : {
                cpp11::sexp res = Rf_allocVector(STRSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, n_samp, IsNamed, IsGmp, lenV,
                               commonLen, STRSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case CPLXSXP : {
                cpp11::sexp res = Rf_allocVector(CPLXSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, n_samp, IsNamed, IsGmp, lenV,
                               commonLen, CPLXSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case RAWSXP : {
                cpp11::sexp res = Rf_allocVector(RAWSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, n_samp, IsNamed, IsGmp, lenV,
                               commonLen, RAWSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case LGLSXP : {
                cpp11::sexp res = Rf_allocVector(LGLSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, n_samp, IsNamed, IsGmp, lenV,
                               commonLen, LGLSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case INTSXP : {
                cpp11::sexp res = Rf_allocVector(INTSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, n_samp, IsNamed, IsGmp, lenV,
                               commonLen, INTSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case REALSXP : {
                cpp11::sexp res = Rf_allocVector(REALSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, n_samp, IsNamed, IsGmp, lenV,
                               commonLen, REALSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } default : {
                cpp11::sexp res = Rf_allocVector(VECSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, n_samp, IsNamed, IsGmp, lenV,
                               commonLen, VECSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            }
        }
    } else {
        cpp11::sexp myList = Rf_allocVector(VECSXP, n_samp);
        SampleApplyFun(myList, v, vectorPass, ptr_vec, mySample,
                       myBigSamp, myReps, stdFun, rho, nthResFun,
                       m, n_samp, IsNamed, IsGmp, lenV);
        return myList;
    }
}

SEXP ApplyFunction(SEXP v, SEXP vectorPass,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   const std::vector<int> &myReps,
                   SEXP stdFun, SEXP rho, SEXP RFunVal,
                   nthResultPtr nthResFun, int m, int n_samp,
                   bool IsNamed, bool IsGmp, int lenV) {

    if (!Rf_isNull(RFunVal)) {
        if (!Rf_isVector(RFunVal)) cpp11::stop("'FUN.VALUE' must be a vector");
        const int commonLen = Rf_length(RFunVal);

        switch (TYPEOF(RFunVal)) {
            case STRSXP : {
                cpp11::sexp res = Rf_allocVector(STRSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, n_samp,
                               IsNamed, IsGmp, lenV, commonLen, STRSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case CPLXSXP : {
                cpp11::sexp res = Rf_allocVector(CPLXSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, n_samp,
                               IsNamed, IsGmp, lenV, commonLen, CPLXSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case RAWSXP : {
                cpp11::sexp res = Rf_allocVector(RAWSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, n_samp,
                               IsNamed, IsGmp, lenV, commonLen, RAWSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case LGLSXP : {
                cpp11::sexp res = Rf_allocVector(LGLSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, n_samp,
                               IsNamed, IsGmp, lenV, commonLen, LGLSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case INTSXP : {
                cpp11::sexp res = Rf_allocVector(INTSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, n_samp,
                               IsNamed, IsGmp, lenV, commonLen, INTSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } case REALSXP : {
                cpp11::sexp res = Rf_allocVector(REALSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, n_samp,
                               IsNamed, IsGmp, lenV, commonLen, REALSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            } default : {
                cpp11::sexp res = Rf_allocVector(VECSXP, n_samp * commonLen);

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, n_samp,
                               IsNamed, IsGmp, lenV, commonLen, VECSXP);

                SetDims(RFunVal, res, commonLen, n_samp);
                return res;
            }
        }
    } else {
        cpp11::sexp myList = Rf_allocVector(VECSXP, n_samp);
        SampleApplyFun(myList, v, vectorPass, mySample, myBigSamp,
                       myReps, stdFun, rho, nthResFun, m, n_samp,
                       IsNamed, IsGmp, lenV);
        return myList;
    }
}

SEXP SampleCombPermApply(SEXP Rv, const std::vector<int> &vInt,
                         const std::vector<double> &vNum,
                         const std::vector<double> &mySample,
                         const std::vector<mpz_class> &myBigSamp,
                         const std::vector<int> &myReps, SEXP stdFun,
                         SEXP rho, SEXP RFunVal, nthResultPtr nthResFun,
                         VecType myType, int n, int m, int sampSize,
                         bool IsNamed, bool IsGmp) {

    switch (myType) {
        case VecType::Character : {
            cpp11::sexp charVec = Rf_duplicate(Rv);
            cpp11::sexp vectorPass = Rf_allocVector(STRSXP, m);

            cpp11::sexp res = ApplyFunction(
                charVec, vectorPass, mySample, myBigSamp, myReps, stdFun,
                rho, RFunVal, nthResFun, m, sampSize, IsNamed, IsGmp, n
            );

            return res;
        } case VecType::Complex : {
            cpp11::sexp vectorPass = Rf_allocVector(CPLXSXP, m);
            Rcomplex* ptr_vec = COMPLEX(vectorPass);
            std::vector<Rcomplex> vCmplx = CppConvert::GetVec<Rcomplex>(Rv);

            cpp11::sexp res = ApplyFunction(
                vCmplx, vectorPass, ptr_vec, mySample, myBigSamp, myReps,
                stdFun, rho, RFunVal, nthResFun, m, sampSize, IsNamed,
                IsGmp, n
            );

            return res;
        } case VecType::Raw : {
            cpp11::sexp vectorPass = Rf_allocVector(RAWSXP, m);
            Rbyte* ptr_vec = RAW(vectorPass);
            std::vector<Rbyte> vByte = CppConvert::GetVec<Rbyte>(Rv);

            cpp11::sexp res = ApplyFunction(
                vByte, vectorPass, ptr_vec, mySample, myBigSamp, myReps,
                stdFun, rho, RFunVal, nthResFun, m, sampSize, IsNamed,
                IsGmp, n
            );

            return res;
        } case VecType::Logical : {
            cpp11::sexp vectorPass = Rf_allocVector(LGLSXP, m);
            int* ptr_vec = LOGICAL(vectorPass);

            cpp11::sexp res = ApplyFunction(
                vInt, vectorPass, ptr_vec, mySample, myBigSamp, myReps,
                stdFun, rho, RFunVal, nthResFun, m, sampSize, IsNamed,
                IsGmp, n
            );

            return res;
        } case VecType::Integer : {
            cpp11::sexp vectorPass = Rf_allocVector(INTSXP, m);
            int* ptr_vec = INTEGER(vectorPass);

            cpp11::sexp res = ApplyFunction(
                vInt, vectorPass, ptr_vec, mySample, myBigSamp, myReps,
                stdFun, rho, RFunVal, nthResFun, m, sampSize, IsNamed,
                IsGmp, n
            );

            return res;
        } default : {
            cpp11::sexp vectorPass = Rf_allocVector(REALSXP, m);
            double* ptr_vec = REAL(vectorPass);

            cpp11::sexp res = ApplyFunction(
                vNum, vectorPass, ptr_vec, mySample, myBigSamp, myReps,
                stdFun, rho, RFunVal, nthResFun, m, sampSize, IsNamed,
                IsGmp, n
            );

            return res;
        }
    }
}
