#include "ClassUtils/PrevCombinatorics.h"
#include "SetUpUtils.h"

void GetPrevious(SEXP mat, SEXP v, std::vector<int> &z,
                 prevIterPtr prevIter, int n, int m, int nRows,
                 const std::vector<int> &freqs, bool IsComb, bool IsMult) {

    const int loc_n1 = IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1);
    const int lastRow = nRows - 1;
    const int loc_m = m;

    // We iterater to the pentultimate row to avoid iterating z one too many times
    for (int count = 0, loc_m1 = m - 1; count < lastRow; ++count) {
        for (int j = 0; j < loc_m; ++j) {
            SET_STRING_ELT(mat, count + j * nRows, STRING_ELT(v, z[j]));
        }

        prevIter(freqs, z, loc_n1, loc_m1);
    }

    // Get the last result
    for (int j = 0; j < loc_m; ++j) {
        SET_STRING_ELT(mat, lastRow + j * nRows, STRING_ELT(v, z[j]));
    }
}

template <typename T>
void GetPrevious(T* mat, const std::vector<T> &v, std::vector<int> &z,
                 prevIterPtr prevIter, int n, int m, int nRows,
                 const std::vector<int> &freqs, bool IsComb, bool IsMult) {

    const int loc_n1 = IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1);
    const int lastRow = nRows - 1;
    const int loc_m = m;

    // We iterater to the pentultimate row to avoid iterating z one too many times
    for (int count = 0, loc_m1 = m - 1; count < lastRow; ++count) {
        for (int j = 0; j < loc_m; ++j) {
            mat[count + j * nRows] = v[z[j]];
        }

        prevIter(freqs, z, loc_n1, loc_m1);
    }

    // Get the last result
    for (int j = 0; j < loc_m; ++j) {
        mat[lastRow + j * nRows] = v[z[j]];
    }
}

SEXP GetPrevCombPerms(SEXP Rv, const std::vector<double> &vNum,
                      const std::vector<int> &vInt,
                      const std::vector<int> &myReps,
                      const std::vector<int> &freqs, std::vector<int> &z,
                      prevIterPtr prevIter, int n, int m, bool IsComb,
                      bool IsMult, int nRows, VecType myType) {

    switch (myType) {
        case VecType::Character : {
            cpp11::sexp charVec = Rf_duplicate(Rv);
            cpp11::sexp res = Rf_allocMatrix(STRSXP, nRows, m);

            GetPrevious(res, charVec, z, prevIter, n, m,
                        nRows, freqs, IsComb, IsMult);

            return res;
        } case VecType::Complex : {
            std::vector<Rcomplex> stlCmplxVec(n);
            Rcomplex* vecCmplx = COMPLEX(Rv);

            for (int i = 0; i < n; ++i) {
                stlCmplxVec[i] = vecCmplx[i];
            }

            cpp11::sexp res = Rf_allocMatrix(CPLXSXP, nRows, m);
            Rcomplex* matCmplx = COMPLEX(res);

            GetPrevious(matCmplx, stlCmplxVec, z, prevIter,
                        n, m, nRows, freqs, IsComb, IsMult);

            return res;
        } case VecType::Raw : {
            std::vector<Rbyte> stlRawVec(n);
            Rbyte* rawVec = RAW(Rv);

            for (int i = 0; i < n; ++i) {
                stlRawVec[i] = rawVec[i];
            }

            cpp11::sexp res = Rf_allocMatrix(RAWSXP, nRows, m);
            Rbyte* rawMat = RAW(res);

            GetPrevious(rawMat, stlRawVec, z, prevIter,
                        n, m, nRows, freqs, IsComb, IsMult);

            return res;
        } case VecType::Logical : {
            std::vector<int> vBool(n, 0);
            int* vecBool = LOGICAL(Rv);

            for (int i = 0; i < n; ++i) {
                vBool[i] = vecBool[i];
            }

            cpp11::sexp res = Rf_allocMatrix(LGLSXP, nRows, m);
            int* matBool = LOGICAL(res);

            GetPrevious(matBool, vBool, z, prevIter, n,
                        m, nRows, freqs, IsComb, IsMult);

            return res;
        } case VecType::Integer : {
            cpp11::sexp res = Rf_allocMatrix(INTSXP, nRows, m);
            int* matInt = INTEGER(res);

            GetPrevious(matInt, vInt, z, prevIter, n,
                        m, nRows, freqs, IsComb, IsMult);

            if (Rf_isFactor(Rv)) {
                SetFactorClass(res, Rv);
            }

            return res;
        } default : {
            cpp11::sexp res = Rf_allocMatrix(REALSXP, nRows, m);
            double* matNum = REAL(res);

            GetPrevious(matNum, vNum, z, prevIter, n,
                        m, nRows, freqs, IsComb, IsMult);

            return res;
        }
    }
}
