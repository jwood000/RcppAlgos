#include "SetUpUtils.h"
#include "FunAssign.h"
#include "NthResult.h"

template <typename T>
void SampleApplyFun(SEXP res, const std::vector<T> &v, SEXP vectorPass,
                    T* ptr_vec, const std::vector<double> &mySample,
                    mpz_t *const myBigSamp, const std::vector<int> &myReps,
                    SEXP func, SEXP rho, nthResultPtr nthResFun, int m,
                    int sampSize, bool IsNamed, bool IsGmp, int lenV,
                    int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);
    SEXP sexpFun = PROTECT(Rf_lang2(func, R_NilValue));

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
        mpz_t mpzDefault;
        mpz_init(mpzDefault);

        for (int count = 0; count < sampSize; ++count) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[count],
                                                 mpzDefault, myReps);

            for (int j = 0; j < m; ++j) {
                ptr_vec[j] = v[z[j]];
            }

            FunAssign(res, vectorPass, sexpFun, rho, commonType,
                      commonLen, count, sampSize, retType);
        }

        mpz_clear(mpzDefault);
    }

    UNPROTECT(1);

    if (IsNamed) {
        SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp);
    }
}

void SampleApplyFun(SEXP res, SEXP v, SEXP vectorPass,
                    const std::vector<double> &mySample,
                    mpz_t *const myBigSamp, const std::vector<int> &myReps,
                    SEXP func, SEXP rho, nthResultPtr nthResFun, int m,
                    int sampSize, bool IsNamed, bool IsGmp, int lenV,
                    int commonLen = 1, int commonType = INTSXP) {

    const int retType = TYPEOF(res);
    SEXP sexpFun = PROTECT(Rf_lang2(func, R_NilValue));

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
        mpz_t mpzDefault;
        mpz_init(mpzDefault);

        for (int count = 0; count < sampSize; ++count) {
            const std::vector<int> z = nthResFun(lenV, m, mySample[count],
                                                 mpzDefault, myReps);

            for (int j = 0; j < m; ++j) {
                SET_STRING_ELT(vectorPass, j, STRING_ELT(v, z[j]));
            }

            FunAssign(res, vectorPass, sexpFun, rho, commonType,
                      commonLen, count, sampSize, retType);
        }

        mpz_clear(mpzDefault);
    }

    UNPROTECT(1);

    if (IsNamed) {
        SetSampleNames(res, IsGmp, sampSize, mySample, myBigSamp);
    }
}

template <typename T>
SEXP ApplyFunction(const std::vector<T> &v, SEXP vectorPass,
                   T* ptr_vec, const std::vector<double> &mySample,
                   mpz_t *const myBigSamp, const std::vector<int> &myReps,
                   SEXP stdFun, SEXP rho, SEXP RFunVal,
                   nthResultPtr nthResFun, int m, int sampSize,
                   bool IsNamed, bool IsGmp, int lenV) {

    if (!Rf_isNull(RFunVal)) {
        if (!Rf_isVector(RFunVal)) cpp11::stop("'FUN.VALUE' must be a vector");
        const int commonLen = Rf_length(RFunVal);

        switch (TYPEOF(RFunVal)) {
            case STRSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(STRSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, sampSize, IsNamed, IsGmp, lenV,
                               commonLen, STRSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case CPLXSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(CPLXSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, sampSize, IsNamed, IsGmp, lenV,
                               commonLen, CPLXSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case RAWSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(RAWSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, sampSize, IsNamed, IsGmp, lenV,
                               commonLen, RAWSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case LGLSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(LGLSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, sampSize, IsNamed, IsGmp, lenV,
                               commonLen, LGLSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case INTSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(INTSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, sampSize, IsNamed, IsGmp, lenV,
                               commonLen, INTSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case REALSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(REALSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, sampSize, IsNamed, IsGmp, lenV,
                               commonLen, REALSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } default : {
                SEXP res = PROTECT(
                    Rf_allocVector(VECSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, ptr_vec, mySample,
                               myBigSamp, myReps, stdFun, rho, nthResFun,
                               m, sampSize, IsNamed, IsGmp, lenV,
                               commonLen, VECSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            }
        }
    } else {
        SEXP myList = PROTECT(Rf_allocVector(VECSXP, sampSize));
        SampleApplyFun(myList, v, vectorPass, ptr_vec, mySample,
                       myBigSamp, myReps, stdFun, rho, nthResFun,
                       m, sampSize, IsNamed, IsGmp, lenV);
        UNPROTECT(1);
        return myList;
    }
}

SEXP ApplyFunction(SEXP v, SEXP vectorPass,
                   const std::vector<double> &mySample,
                   mpz_t *const myBigSamp, const std::vector<int> &myReps,
                   SEXP stdFun, SEXP rho, SEXP RFunVal,
                   nthResultPtr nthResFun, int m, int sampSize,
                   bool IsNamed, bool IsGmp, int lenV) {

    if (!Rf_isNull(RFunVal)) {
        if (!Rf_isVector(RFunVal)) cpp11::stop("'FUN.VALUE' must be a vector");
        const int commonLen = Rf_length(RFunVal);

        switch (TYPEOF(RFunVal)) {
            case STRSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(STRSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, sampSize,
                               IsNamed, IsGmp, lenV, commonLen, STRSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case CPLXSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(CPLXSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, sampSize,
                               IsNamed, IsGmp, lenV, commonLen, CPLXSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case RAWSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(RAWSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, sampSize,
                               IsNamed, IsGmp, lenV, commonLen, RAWSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case LGLSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(LGLSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, sampSize,
                               IsNamed, IsGmp, lenV, commonLen, LGLSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case INTSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(INTSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, sampSize,
                               IsNamed, IsGmp, lenV, commonLen, INTSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } case REALSXP : {
                SEXP res = PROTECT(
                    Rf_allocVector(REALSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, sampSize,
                               IsNamed, IsGmp, lenV, commonLen, REALSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            } default : {
                SEXP res = PROTECT(
                    Rf_allocVector(VECSXP, sampSize * commonLen)
                );

                SampleApplyFun(res, v, vectorPass, mySample, myBigSamp,
                               myReps, stdFun, rho, nthResFun, m, sampSize,
                               IsNamed, IsGmp, lenV, commonLen, VECSXP);

                SetDims(RFunVal, res, commonLen, sampSize);
                UNPROTECT(1);
                return res;
            }
        }
    } else {
        SEXP myList = PROTECT(Rf_allocVector(VECSXP, sampSize));
        SampleApplyFun(myList, v, vectorPass, mySample, myBigSamp,
                       myReps, stdFun, rho, nthResFun, m, sampSize,
                       IsNamed, IsGmp, lenV);
        UNPROTECT(1);
        return myList;
    }
}

SEXP SampleCombPermApply(SEXP Rv, const std::vector<int> &vInt,
                         const std::vector<double> &vNum,
                         const std::vector<double> &mySample,
                         mpz_t *const myBigSamp,
                         const std::vector<int> &myReps, SEXP stdFun,
                         SEXP rho, SEXP RFunVal, nthResultPtr nthResFun,
                         VecType myType, int n, int m, int sampSize,
                         bool IsNamed, bool IsGmp) {

    switch (myType) {
        case VecType::Character : {
            SEXP charVec = PROTECT(Rf_duplicate(Rv));
            SEXP vectorPass = PROTECT(Rf_allocVector(STRSXP, m));

            SEXP res = PROTECT(
                ApplyFunction(charVec, vectorPass, mySample, myBigSamp,
                              myReps, stdFun, rho, RFunVal, nthResFun,
                              m, sampSize, IsNamed, IsGmp, n)
            );

            MpzClearVec(myBigSamp, sampSize, IsGmp);
            UNPROTECT(3);
            return res;
        } case VecType::Complex : {
            SEXP vectorPass = PROTECT(Rf_allocVector(CPLXSXP, m));
            Rcomplex* ptr_vec = COMPLEX(vectorPass);

            Rcomplex* cmplxVec = COMPLEX(Rv);
            std::vector<Rcomplex> vCmplx(cmplxVec, cmplxVec + n);

            SEXP res = PROTECT(
                ApplyFunction(vCmplx, vectorPass, ptr_vec, mySample,
                              myBigSamp, myReps, stdFun, rho, RFunVal,
                              nthResFun, m, sampSize, IsNamed, IsGmp, n)
            );

            MpzClearVec(myBigSamp, sampSize, IsGmp);
            UNPROTECT(2);
            return res;
        } case VecType::Raw : {
            SEXP vectorPass = PROTECT(Rf_allocVector(RAWSXP, m));
            Rbyte* ptr_vec = RAW(vectorPass);

            Rbyte* rawVec = RAW(Rv);
            std::vector<Rbyte> vByte(rawVec, rawVec + n);

            SEXP res = PROTECT(
                ApplyFunction(vByte, vectorPass, ptr_vec, mySample,
                              myBigSamp, myReps, stdFun, rho, RFunVal,
                              nthResFun, m, sampSize, IsNamed, IsGmp, n)
            );

            MpzClearVec(myBigSamp, sampSize, IsGmp);
            UNPROTECT(2);
            return res;
        } case VecType::Logical : {
            SEXP vectorPass = PROTECT(Rf_allocVector(LGLSXP, m));
            int* ptr_vec = LOGICAL(vectorPass);

            SEXP res = PROTECT(
                ApplyFunction(vInt, vectorPass, ptr_vec, mySample,
                              myBigSamp, myReps, stdFun, rho, RFunVal,
                              nthResFun, m, sampSize, IsNamed, IsGmp, n)
            );

            MpzClearVec(myBigSamp, sampSize, IsGmp);
            UNPROTECT(2);
            return res;
        } case VecType::Integer : {
            SEXP vectorPass = PROTECT(Rf_allocVector(INTSXP, m));
            int* ptr_vec = INTEGER(vectorPass);

            SEXP res = PROTECT(
                ApplyFunction(vInt, vectorPass, ptr_vec, mySample,
                              myBigSamp, myReps, stdFun, rho, RFunVal,
                              nthResFun, m, sampSize, IsNamed, IsGmp, n)
            );

            MpzClearVec(myBigSamp, sampSize, IsGmp);
            UNPROTECT(2);
            return res;
        } default : {
            SEXP vectorPass = PROTECT(Rf_allocVector(REALSXP, m));
            double* ptr_vec = REAL(vectorPass);

            SEXP res = PROTECT(
                ApplyFunction(vNum, vectorPass, ptr_vec, mySample,
                              myBigSamp, myReps, stdFun, rho, RFunVal,
                              nthResFun, m, sampSize, IsNamed, IsGmp, n)
            );

            MpzClearVec(myBigSamp, sampSize, IsGmp);
            UNPROTECT(2);
            return res;
        }
    }
}
