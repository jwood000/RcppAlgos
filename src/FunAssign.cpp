#include "FunAssign.h"
#include "cpp11/protect.hpp"
#include "cpp11/sexp.hpp"

// See function do_vapply found here:
//     https://github.com/wch/r-source/blob/trunk/src/main/apply.c
void VapplyAssign(SEXP ans, SEXP vectorPass,
                  SEXP sexpFun, SEXP rho, int commonType,
                  int commonLen, int count, int nRows) {

    SEXPTYPE valType;
    PROTECT_INDEX indx;
    SETCADR(sexpFun, vectorPass);

    SEXP val = Rf_eval(sexpFun, rho);
    PROTECT_WITH_INDEX(val, &indx);

    if (Rf_length(val) != commonLen) {
        cpp11::stop("values must be length %d,\n but "
                 "FUN(X[[%d]]) result is length %d",
                 commonLen, count + 1, Rf_length(val));
    }

    valType = TYPEOF(val);

    if (static_cast<int>(valType) != commonType) {
        bool okay = false;
        switch (commonType) {
            case CPLXSXP: okay = (valType == REALSXP) || (valType == INTSXP)
                || (valType == LGLSXP); break;
            case REALSXP: okay = (valType == INTSXP) || (valType == LGLSXP); break;
            case INTSXP:  okay = (valType == LGLSXP); break;
        }

        if (!okay) {
            cpp11::stop("values must be type '%s',\n but FUN(X[[%d]]) result is type '%s'",
                  Rf_type2char(commonType), count + 1, Rf_type2char(valType));
        }

        REPROTECT(val = Rf_coerceVector(val, commonType), indx);
    }

    if(commonLen == 1) { // common case
        switch (commonType) {
            case CPLXSXP: COMPLEX(ans)[count] = COMPLEX(val)[0]; break;
            case REALSXP: REAL(ans)   [count] = REAL   (val)[0]; break;
            case INTSXP:  INTEGER(ans)[count] = INTEGER(val)[0]; break;
            case LGLSXP:  LOGICAL(ans)[count] = LOGICAL(val)[0]; break;
            case RAWSXP:  RAW(ans)    [count] = RAW    (val)[0]; break;
            case STRSXP:  SET_STRING_ELT(ans, count, STRING_ELT(val, 0)); break;
            case VECSXP:  SET_VECTOR_ELT(ans, count, VECTOR_ELT(val, 0)); break;
        }
    } else { // commonLen > 1 (typically, or == 0) :

        const std::size_t unRows = nRows;
        const std::size_t uCommonLen = commonLen;
        const std::size_t uCount = count;

        switch (commonType) {
            case REALSXP: {
                double* ans_dbl = REAL(ans);
                double* val_dbl = REAL(val);

                for (std::size_t j = 0; j < uCommonLen; j++) {
                    ans_dbl[uCount + j * unRows] = val_dbl[j];
                }

                break;
            } case INTSXP: {
                int* ans_int = INTEGER(ans);
                int* val_int = INTEGER(val);

                for (std::size_t j = 0; j < uCommonLen; j++) {
                    ans_int[uCount + j * unRows] = val_int[j];
                }

                break;
            } case LGLSXP: {
                int* ans_bool = LOGICAL(ans);
                int* val_bool = LOGICAL(val);

                for (std::size_t j = 0; j < uCommonLen; j++) {
                    ans_bool[uCount + j * unRows] = val_bool[j];
                }

                break;
            } case RAWSXP: {
                Rbyte* ans_raw = RAW(ans);
                Rbyte* val_raw = RAW(val);

                for (std::size_t j = 0; j < uCommonLen; j++) {
                    ans_raw[uCount + j * unRows] = val_raw[j];
                }

                break;
            } case CPLXSXP: {
                Rcomplex* ans_cmplx = COMPLEX(ans);
                Rcomplex* val_cmplx = COMPLEX(val);

                for (std::size_t j = 0; j < uCommonLen; j++) {
                    ans_cmplx[uCount + j * unRows] = val_cmplx[j];
                }

                break;
            } case STRSXP: {
                for (std::size_t j = 0; j < uCommonLen; j++) {
                    SET_STRING_ELT(ans, uCount + j * unRows,
                                   STRING_ELT(val, j));
                }

                break;
            } case VECSXP: {
                for (std::size_t j = 0; j < uCommonLen; j++) {
                    SET_VECTOR_ELT(ans, uCount + j * unRows,
                                   VECTOR_ELT(val, j));
                }

                break;
            }
        }
    }

    UNPROTECT(1);
}

void FunAssign(SEXP res, SEXP vectorPass, SEXP sexpFun,
               SEXP rho, int commonType, int commonLen,
               int count, int nRows, int retType) {

    if (retType == VECSXP) {
        SETCADR(sexpFun, Rf_duplicate(vectorPass));
        SET_VECTOR_ELT(res, count, Rf_eval(sexpFun, rho));
    } else {
        VapplyAssign(res, vectorPass, sexpFun, rho,
                     commonType, commonLen, count, nRows);
    }
}

void SetDims(SEXP RFunVal, SEXP res, int commonLen, int nRows) {

    cpp11::sexp dim_v = Rf_getAttrib(RFunVal, R_DimSymbol);
    const bool array_value = (TYPEOF(dim_v) == INTSXP && LENGTH(dim_v) >= 1);

    if (commonLen != 1) {
        const int rnk_v = array_value ? LENGTH(dim_v) : 1;
        cpp11::sexp dim = Rf_allocVector(INTSXP, rnk_v + 1);
        INTEGER(dim)[0] = nRows;

        if(array_value) {
            for(int j = 0; j < rnk_v; j++) {
                INTEGER(dim)[j + 1] = INTEGER(dim_v)[j];
            }
        } else {
            INTEGER(dim)[rnk_v] = commonLen;
        }

        Rf_setAttrib(res, R_DimSymbol, dim);
    }

}

