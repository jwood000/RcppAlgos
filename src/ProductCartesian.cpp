#include "cpp11/logicals.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/strings.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/list.hpp"

#include "SetUpUtils.h"
#include <algorithm>
#include <cstdint>
#include <numeric>

double CartesianCount(const std::vector<int> &lenGrps) {
    return std::accumulate(lenGrps.begin(), lenGrps.end(),
                           1.0, std::multiplies<double>());
}

void CartesianCountGmp(mpz_class &result, const std::vector<int> &lenGrps) {

    result = 1;

    for (auto len: lenGrps) {
        result *= len;
    }
}

std::vector<int> nthProduct(double dblIdx, const std::vector<int> &lenGrp) {

    double index1 = dblIdx;
    const int m = lenGrp.size();

    std::vector<int> res(m);
    double temp = CartesianCount(lenGrp);

    for (int k = 0; k < m; ++k) {
        temp /= lenGrp[k];
        int j = static_cast<int>(index1 / temp);
        res[k] = j;
        index1 -= (temp * j);
    }

    return res;
}

bool nextProduct(const std::vector<int> &lenGrps,
                 std::vector<int> &z, int m) {

    if (z.back() < (lenGrps.back() - 1)) {
        ++z.back();
        return true;
    } else {
        z.back() = 0;

        for (int i = m - 2; i >= 0; --i) {
            if (z[i] < (lenGrps[i] - 1)) {
                ++z[i];
                return true;
            } else {
                z[i] = 0;
            }
        }
    }

    return false;
}

template <typename T>
void GetPureOutput(T* mat,
                   const std::vector<int> &mat_idx,
                   const std::vector<int> &lenGrps,
                   const std::vector<T> &standardVec,
                   std::vector<int> z, int nCols, int nRows) {

    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            mat[i + j * nRows] = standardVec[mat_idx[j + z[j] * nCols]];
        }

        nextProduct(lenGrps, z, nCols);
    }
}

void GetCharOutput(cpp11::writable::strings_matrix<> &mat,
                   const std::vector<std::vector<int>> &idx,
                   const cpp11::strings &charVec,
                   int nCols, int nRows) {

    const int lastCol = nCols - 1;

    for (int i = 0, preSecLen = nRows; i < lastCol; ++i) {
        const int grpSize = idx[i].size();
        const int secLen = preSecLen / grpSize;
        int strt = 0;

        for (auto v_i: idx[i]) {
            SEXP val = PROTECT(STRING_ELT(charVec, v_i));

            for (int j = 0; j < nRows; j += preSecLen) {
                for (int k = 0, q = strt, m_idx = i * nRows + j;
                     k < secLen; ++k, ++q) {
                    SET_STRING_ELT(mat, m_idx + q, val);
                }
            }

            UNPROTECT(1);
            strt += secLen;
        }

        preSecLen /= grpSize;
    }

    const int lastGrpSize = idx.back().size();
    int strt = 0;

    for (auto v_i: idx.back()) {
        SEXP val = PROTECT(STRING_ELT(charVec, v_i));

        for (int k = strt, m_idx = lastCol * nRows;
             k < nRows; k += lastGrpSize) {
            SET_STRING_ELT(mat, m_idx + k, val);
        }

        UNPROTECT(1);
        ++strt;
    }
}

template <typename T>
void GetPureOutput(T* mat, const std::vector<std::vector<int>> &idx,
                   const std::vector<T> &standardVec, int nCols, int nRows) {

    const int lastCol = nCols - 1;

    for (int i = 0, preSecLen = nRows; i < lastCol; ++i) {
        const int grpSize = idx[i].size();
        const int secLen = preSecLen / grpSize;
        int strt = 0;

        for (auto v_i: idx[i]) {
            auto&& val = standardVec[v_i];

            for (int j = 0; j < nRows; j += preSecLen) {
                for (int k = 0, q = strt, m_idx = i * nRows + j;
                     k < secLen; ++k, ++q) {
                    mat[m_idx + q] = val;
                }
            }

            strt += secLen;
        }

        preSecLen /= grpSize;
    }

    const int lastGrpSize = idx.back().size();
    int strt = 0;

    for (auto v_i: idx.back()) {
        auto&& val = standardVec[v_i];

        for (int k = strt, m_idx = lastCol * nRows;
             k < nRows; k += lastGrpSize) {
            mat[m_idx + k] = val;
        }

        ++strt;
    }
}

template <typename T>
void PoulateColumn(T* vec,
                   const std::vector<int> &idx,
                   const std::vector<T> &poolVec,
                   int nCols, int nRows, int preSecLen, int i) {

    const int grpSize = idx.size();
    int strt = 0;

    if (i < (nCols - 1)) {
        const int secLen = preSecLen / grpSize;

        for (auto v_i: idx) {
            auto&& val = poolVec[v_i];

            for (int j = 0; j < nRows; j += preSecLen) {
                for (int k = 0, q = strt; k < secLen; ++k, ++q) {
                    vec[j + q] = val;
                }
            }

            strt += secLen;
        }
    } else {
        for (auto v_i: idx) {
            auto&& val = poolVec[v_i];

            for (int k = strt; k < nRows; k += grpSize) {
                vec[k] = val;
            }

            ++strt;
        }
    }
}


template <typename T>
void GroupsMain(T* GroupsMat, SEXP res, nextGrpFunc nextCmbGrp,
                nthFuncDbl nthCmbGrp, nthFuncGmp nthCmbGrpGmp,
                finalTouchFunc FinalTouch, const std::vector<T> &v,
                std::vector<int> z, const std::vector<double> &mySample,
                const std::vector<mpz_class> &myBigSamp, mpz_class lowerMpz,
                double lower, int n, int nRows, int nThreads, bool IsArray,
                bool IsNamed, bool Parallel, bool IsSample, bool IsGmp) {

    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(GroupsMat, nRows, n);
        std::vector<std::thread> threads;

        int step = 0;
        int stepSize = nRows / nThreads;
        int nextStep = stepSize;

        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize,
             nextStep += stepSize) {

            threads.emplace_back(
                std::cref(ParallelGlue<T>), std::ref(parMat), std::cref(v),
                nextCmbGrp, nthCmbGrp, nthCmbGrpGmp, std::cref(mySample),
                std::cref(myBigSamp), z, n, step, nextStep, IsSample, IsGmp
            );

            GetStartGrp(nthCmbGrp, nthCmbGrpGmp, z,
                        lowerMpz, lower, stepSize, IsGmp);
        }

        threads.emplace_back(
            std::cref(ParallelGlue<T>), std::ref(parMat), std::cref(v),
            nextCmbGrp, nthCmbGrp, nthCmbGrpGmp, std::cref(mySample),
            std::cref(myBigSamp), z, n, step, nRows, IsSample, IsGmp
        );

        for (auto& thr: threads) {
            thr.join();
        }

        FinalTouch(res, IsArray, nRows, IsNamed,
                   mySample, myBigSamp, IsSample);
    } else {
        SerialGlue(GroupsMat, res, v, nextCmbGrp, nthCmbGrp, nthCmbGrpGmp,
                   FinalTouch, mySample, myBigSamp, z, n, nRows, IsArray,
                   IsSample, IsNamed, IsGmp);
    }
}

SEXP GlueProdCart(
    const std::vector<std::vector<int>> &idx, const std::vector<int> &mat_idx,
    const std::vector<int> &typeCheck, const std::vector<int> &IsFactor,
    const cpp11::list &RList, const std::vector<int> &intVec,
    const std::vector<double> &dblVec, const std::vector<int> &boolVec,
    const std::vector<Rcomplex> &cmplxVec, const std::vector<Rbyte> &rawVec,
    const cpp11::strings &charVec, const std::vector<int> &lenGrps,
    std::vector<int> &z, int nRows, int nCols, bool IsDF
) {

    if (IsDF) {
        cpp11::writable::list DataFrame(nCols);

        for (int i = 0, preSecLen = nRows; i < nCols; ++i) {
            switch (TYPEOF(RList[i])) {
                case INTSXP: {
                    cpp11::sexp res = Rf_allocVector(INTSXP, nRows);
                    int* intSexpVec = INTEGER(res);

                    if (IsFactor[i]) {
                        PoulateColumn(intSexpVec, idx[i], intVec,
                                      nCols, nRows, preSecLen, i);
                        SetFactorClass(res, RList[i]);
                    } else {
                        PoulateColumn(intSexpVec, idx[i], intVec,
                                      nCols, nRows, preSecLen, i);
                    }

                    DataFrame[i] = res;
                    break;
                } case LGLSXP: {
                    cpp11::sexp res = Rf_allocVector(LGLSXP, nRows);
                    int* boolSexpVec = LOGICAL(res);
                    PoulateColumn(boolSexpVec, idx[i], boolVec,
                                  nCols, nRows, preSecLen, i);
                    DataFrame[i] = res;
                    break;
                } case REALSXP: {
                    cpp11::sexp res = Rf_allocVector(REALSXP, nRows);
                    double* dblSexpVec = REAL(res);
                    PoulateColumn(dblSexpVec, idx[i], dblVec,
                                  nCols, nRows, preSecLen, i);
                    DataFrame[i] = res;
                    break;
                } case CPLXSXP : {
                    cpp11::sexp res = Rf_allocVector(CPLXSXP, nRows);
                    Rcomplex* cmplxSexpVec = COMPLEX(res);
                    PoulateColumn(cmplxSexpVec, idx[i], cmplxVec,
                                  nCols, nRows, preSecLen, i);
                    DataFrame[i] = res;
                    break;
                } case RAWSXP : {
                    cpp11::sexp res = Rf_allocVector(RAWSXP, nRows);
                    Rbyte* rawSexpVec = RAW(res);
                    PoulateColumn(rawSexpVec, idx[i], rawVec,
                                  nCols, nRows, preSecLen, i);
                    DataFrame[i] = res;
                    break;
                } case STRSXP: {
                    cpp11::writable::strings sexpVec(nRows);
                    const int grpSize = idx[i].size();
                    int strt = 0;

                    if (i < (nCols - 1)) {
                        const int secLen = preSecLen / grpSize;

                        for (auto v_i: idx[i]) {
                            SEXP val = PROTECT(STRING_ELT(charVec, v_i));

                            for (int j = 0; j < nRows; j += preSecLen) {
                                for (int k = 0, q = strt;
                                     k < secLen; ++k, ++q) {
                                    SET_STRING_ELT(sexpVec, j + q, val);
                                }
                            }

                            UNPROTECT(1);
                            strt += secLen;
                        }
                    } else {
                        for (auto v_i: idx[i]) {
                            SEXP val = PROTECT(STRING_ELT(charVec, v_i));

                            for (int k = strt; k < nRows; k += grpSize) {
                                SET_STRING_ELT(sexpVec, k, val);
                            }

                            UNPROTECT(1);
                            ++strt;
                        }
                    }

                    DataFrame[i] = sexpVec;
                    break;
                }
            }

            preSecLen /= idx[i].size();
        }

        DataFrame.attr("row.names") = {NA_INTEGER, -nRows};
        DataFrame.names() = RList.names();
        DataFrame.attr("class") = "data.frame";
        return DataFrame;
    } else {
        switch (TYPEOF(RList[0])) {
            case INTSXP : {
                cpp11::sexp res = Rf_allocMatrix(INTSXP, nRows, nCols);
                int* intMat = INTEGER(res);
                // GetPureOutput(intMat, idx, intVec, nCols, nRows);
                GetPureOutput(intMat, mat_idx, lenGrps, intVec, z, nCols, nRows);
                if (typeCheck[tFac]) SetFactorClass(res, RList[0]);
                return res;
            } case LGLSXP : {
                cpp11::sexp res = Rf_allocMatrix(LGLSXP, nRows, nCols);
                int* boolMat = LOGICAL(res);
                GetPureOutput(boolMat, idx, boolVec, nCols, nRows);
                return res;
            } case RAWSXP : {
                cpp11::sexp res = Rf_allocMatrix(RAWSXP, nRows, nCols);
                Rbyte* rawMat = RAW(res);
                GetPureOutput(rawMat, idx, rawVec, nCols, nRows);
                return res;
            } case CPLXSXP : {
                cpp11::sexp res = Rf_allocMatrix(CPLXSXP, nRows, nCols);
                Rcomplex* cmplxMat = COMPLEX(res);
                GetPureOutput(cmplxMat, idx, cmplxVec, nCols, nRows);
                return res;
            } case REALSXP : {
                cpp11::sexp res = Rf_allocMatrix(REALSXP, nRows, nCols);
                double* dblMat = REAL(res);
                GetPureOutput(dblMat, idx, dblVec, nCols, nRows);
                return res;
            } case STRSXP : {
                cpp11::writable::strings_matrix<> charMat(nRows, nCols);
                GetCharOutput(charMat, idx, charVec, nCols, nRows);
                return charMat;
            }
        }
    }
}

[[cpp11::register]]
SEXP ExpandGridCpp(cpp11::list RList, SEXP Rlow, SEXP Rhigh,
                   SEXP RNumThreads, SEXP RmaxThreads) {

    const int nCols = Rf_length(RList);
    std::vector<int> IsFactor(nCols);
    std::vector<int> lenGrps(nCols);

    for (int i = 0; i < nCols; ++i) {
        if (Rf_isFactor(RList[i])) {
            IsFactor[i] = 1;
        } else {
            IsFactor[i] = 0;
        }

        lenGrps[i] = Rf_length(RList[i]);
    }

    const int sumLength = std::accumulate(
        lenGrps.begin(), lenGrps.end(), 0
    );

    std::vector<std::vector<int>> myVec(nCols);
    std::vector<int> typeCheck(N_TYPES, 0);

    cpp11::writable::strings charVec(sumLength);
    std::vector<Rcomplex> cmplxVec(sumLength);
    std::vector<Rbyte> rawVec(sumLength);
    std::vector<double> dblVec(sumLength);
    std::vector<int> intVec(sumLength);
    std::vector<int> boolVec(sumLength);

    VecType myType = VecType::Integer;

    for (int i = 0, strt = 0; i < nCols; ++i) {
        switch(TYPEOF(RList[i])) {
            case INTSXP : {
                if (IsFactor[i]) {
                    typeCheck[tFac] = 1;
                } else {
                    typeCheck[tInt] = 1;
                }

                std::vector<int> temp = CppConvert::GetVec<int>(RList[i]);
                std::copy(temp.begin(), temp.end(), intVec.begin() + strt);
                myType = VecType::Integer;
                break;
            } case LGLSXP : {
                std::vector<int> temp = CppConvert::GetVec<int>(RList[i]);
                std::copy(temp.begin(), temp.end(), boolVec.begin() + strt);
                typeCheck[tLog] = 1;
                myType = VecType::Logical;
                break;
            } case CPLXSXP : {
                std::vector<Rcomplex> temp =
                    CppConvert::GetVec<Rcomplex>(RList[i]);
                std::copy(temp.begin(), temp.end(), cmplxVec.begin() + strt);
                typeCheck[tCpx] = 1;
                myType = VecType::Complex;
                break;
            } case RAWSXP : {
                std::vector<Rbyte> temp = CppConvert::GetVec<Rbyte>(RList[i]);
                std::copy(temp.begin(), temp.end(), rawVec.begin() + strt);
                typeCheck[tRaw] = 1;
                myType = VecType::Raw;
                break;
            } case REALSXP : {
                std::vector<double> temp =
                    CppConvert::GetVec<double>(RList[i]);
                std::copy(temp.begin(), temp.end(), dblVec.begin() + strt);
                typeCheck[tDbl] = 1;
                myType = VecType::Numeric;
                break;
            } case STRSXP : {
                for (int j = 0; j < lenGrps[i]; ++j) {
                    charVec[strt + j] = STRING_ELT(RList[i], j);
                }

                typeCheck[tStr] = 1;
                myType = VecType::Character;
                break;
            }
        }

        std::vector<int> idx(lenGrps[i]);
        std::iota(idx.begin(), idx.end(), strt);

        myVec[i] = idx;
        strt += lenGrps[i];
    }

    int mySum = std::accumulate(typeCheck.cbegin(), typeCheck.cend(), 0);

    // We need to check to see if there is overlap in factor levels
    if (typeCheck[tFac] && mySum == 1) {
        std::vector<std::string> testLevels;

        for (int i = 0; i < nCols; ++i) {
            if (IsFactor[i]) {
                cpp11::strings facVec(Rf_getAttrib(RList[i], R_LevelsSymbol));
                std::vector<std::string> strVec;

                for (auto f: facVec) {
                    const std::string temp(CHAR(f));
                    strVec.push_back(temp);
                }

                if (testLevels.size()) {
                    if (strVec != testLevels) {
                        ++mySum;
                        break;
                    }
                } else {
                    testLevels = strVec;
                }
            }
        }
    }

    bool IsDF = (mySum > 1) ? true : false;
    int nRows = 0;
    int nThreads = 1;
    int maxThreads = 1;

    const double computedRows = CartesianCount(lenGrps);
    const bool IsGmp = (computedRows > Significand53);

    mpz_class computedRowsMpz;

    if (IsGmp) {
        CartesianCountGmp(computedRowsMpz, lenGrps);
    }

    double lower = 0;
    double upper = 0;
    bool Parallel = false;

    bool bLower = false;
    bool bUpper = false;

    mpz_class lowerMpz;
    mpz_class upperMpz;

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz, upperMpz, computedRowsMpz, computedRows);

    double userNumRows = 0;
    SetNumResults(IsGmp, bLower, bUpper, true, upperMpz,
                  lowerMpz, lower, upper, computedRows,
                  computedRowsMpz, nRows, userNumRows);

    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows,
               myType, nThreads, RNumThreads, limit);
    std::vector<int> startZ = nthProduct(lower, lenGrps);

    const int maxLen = *std::max_element(lenGrps.begin(), lenGrps.end());
    std::vector<int> mat_idx(maxLen * nCols);

    // transposing myVec so that each row represents
    // each element of the given list
    for (int i = 0; i < nCols; ++i) {
        for (int j = 0; j < lenGrps[i]; ++j) {
            mat_idx[i + j * nCols] = myVec[i][j];
        }
    }

    cpp11::sexp res = GlueProdCart(
        myVec, mat_idx, typeCheck, IsFactor, RList, intVec, dblVec, boolVec,
        cmplxVec, rawVec, charVec, lenGrps, startZ, nRows, nCols, IsDF
    );

    return res;
}