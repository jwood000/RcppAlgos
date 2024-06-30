#include "cpp11/logicals.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/doubles.hpp"
#include "cpp11/strings.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/list.hpp"

#include "SetUpUtils.h"
#include "RMatrix.h"
#include <functional>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <thread>

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

    for (auto &v_i: res) {
        v_i *= m;
    }

    return res;
}

bool nextProduct(const std::vector<int> &lenGrps,
                 std::vector<int> &z, int m) {

    if (z.back() < lenGrps.back()) {
        z.back() += m;
        return true;
    } else {
        z.back() = 0;

        for (int i = m - 2; i >= 0; --i) {
            if (z[i] < lenGrps[i]) {
                z[i] += m;
                return true;
            } else {
                z[i] = 0;
            }
        }
    }

    return false;
}

template <typename T>
void GetPureOutput(T* mat, const std::vector<int> &idx,
                   const std::vector<int> &lenGrps,
                   const std::vector<T> &standardVec,
                   std::vector<int> &z, int nCols, int nRows) {

    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            mat[i + j * nRows] = standardVec[idx[j + z[j]]];
        }

        nextProduct(lenGrps, z, nCols);
    }
}

template <typename T>
void ParallelProduct(
    RcppParallel::RMatrix<T> &mat, const std::vector<int> &idx,
    const std::vector<int> &lenGrps, const std::vector<T> &standardVec,
    std::vector<int> z, int nCols, int strt, int nRows
) {

    for (int i = strt; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            mat(i, j) = standardVec[idx[j + z[j]]];
        }

        nextProduct(lenGrps, z, nCols);
    }
}

void GetCharOutput(cpp11::writable::strings_matrix<> &mat,
                   const std::vector<int> &idx,
                   const std::vector<int> &lenGrps,
                   const cpp11::strings &charVec,
                   std::vector<int> &z, int nCols, int nRows) {

    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            SET_STRING_ELT(mat, i + j * nRows,
                           STRING_ELT(charVec, idx[j + z[j]]));
        }

        nextProduct(lenGrps, z, nCols);
    }
}

template <typename T>
void PopulateColumn(T* vec,
                    const std::vector<int> &idx,
                    const std::vector<int> &all_idx,
                    const std::vector<T> &poolVec,
                    int nCols, int nRows, int col_idx) {

    for (int i = 0; i < nRows; ++i) {
        vec[i] = poolVec[idx[col_idx + all_idx[i * nCols + col_idx]]];
    }
}

template <typename T>
void PureOutputMain(
    T* mat, const std::vector<int> &idx,
    const std::vector<int> &lenGrps, const std::vector<T> &standardVec,
    std::vector<int> &z, int nCols, int nRows, int nThreads, bool Parallel
) {

    if (Parallel) {
        std::vector<int> lenNxtPr(lenGrps);

        for (auto &v_i: lenNxtPr) {
            v_i = (v_i / nCols) + 1;
        }

        RcppParallel::RMatrix<T> parMat(mat, nRows, nCols);
        std::vector<std::thread> threads;

        int step = 0;
        int stepSize = nRows / nThreads;
        int nextStep = stepSize;

        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize,
             nextStep += stepSize) {

            threads.emplace_back(
                std::cref(ParallelProduct<T>), std::ref(parMat),
                std::cref(idx), std::cref(lenGrps), std::cref(standardVec),
                z, nCols, step, nextStep
            );

            z = nthProduct(nextStep, lenNxtPr);
        }

        threads.emplace_back(
            std::cref(ParallelProduct<T>), std::ref(parMat),
            std::cref(idx), std::cref(lenGrps), std::cref(standardVec),
            z, nCols, step, nRows
        );

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        GetPureOutput(mat, idx, lenGrps, standardVec, z, nCols, nRows);
    }
}

SEXP GlueProdCart(
    const std::vector<int> &idx, const std::vector<int> &typeCheck,
    const std::vector<int> &IsFactor, const cpp11::list &RList,
    const std::vector<int> &intVec, const std::vector<double> &dblVec,
    const std::vector<int> &boolVec, const std::vector<Rcomplex> &cmplxVec,
    const std::vector<Rbyte> &rawVec, const cpp11::strings &charVec,
    const std::vector<int> &lenGrps, std::vector<int> &z,
    int nRows, int nCols, bool IsDF, int nThreads, bool Parallel
) {

    if (IsDF) {
        cpp11::writable::list DataFrame(nCols);
        std::vector<int> all_idx(nRows * nCols);

        for (int i = 0; i < nRows; ++i) {
            std::copy(z.begin(), z.end(), all_idx.begin() + i * nCols);
            nextProduct(lenGrps, z, nCols);
        }

        for (int j = 0; j < nCols; ++j) {
            switch (TYPEOF(RList[j])) {
                case INTSXP: {
                    cpp11::sexp res = Rf_allocVector(INTSXP, nRows);
                    int* intSexpVec = INTEGER(res);

                    PopulateColumn(intSexpVec, idx, all_idx,
                                   intVec, nCols, nRows, j);

                    if (IsFactor[j]) SetFactorClass(res, RList[j]);
                    DataFrame[j] = res;
                    break;
                } case LGLSXP: {
                    cpp11::sexp res = Rf_allocVector(LGLSXP, nRows);
                    int* boolSexpVec = LOGICAL(res);

                    PopulateColumn(boolSexpVec, idx, all_idx,
                                   boolVec, nCols, nRows, j);

                    DataFrame[j] = res;
                    break;
                } case REALSXP: {
                    cpp11::sexp res = Rf_allocVector(REALSXP, nRows);
                    double* dblSexpVec = REAL(res);

                    PopulateColumn(dblSexpVec, idx, all_idx,
                                   dblVec, nCols, nRows, j);

                    DataFrame[j] = res;
                    break;
                } case CPLXSXP : {
                    cpp11::sexp res = Rf_allocVector(CPLXSXP, nRows);
                    Rcomplex* cmplxSexpVec = COMPLEX(res);

                    PopulateColumn(cmplxSexpVec, idx, all_idx,
                                   cmplxVec, nCols, nRows, j);

                    DataFrame[j] = res;
                    break;
                } case RAWSXP : {
                    cpp11::sexp res = Rf_allocVector(RAWSXP, nRows);
                    Rbyte* rawSexpVec = RAW(res);

                    PopulateColumn(rawSexpVec, idx, all_idx,
                                   rawVec, nCols, nRows, j);

                    DataFrame[j] = res;
                    break;
                } case STRSXP: {
                    cpp11::writable::strings sexpVec(nRows);

                    for (int i = 0; i < nRows; ++i) {
                        SET_STRING_ELT(
                            sexpVec, i,
                            STRING_ELT(charVec, idx[j + all_idx[i * nCols + j]])
                        );
                    }

                    DataFrame[j] = sexpVec;
                    break;
                }
            }
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
                PureOutputMain(intMat, idx, lenGrps, intVec,
                               z, nCols, nRows, nThreads, Parallel);
                if (typeCheck[tFac]) SetFactorClass(res, RList[0]);
                return res;
            } case LGLSXP : {
                cpp11::sexp res = Rf_allocMatrix(LGLSXP, nRows, nCols);
                int* boolMat = LOGICAL(res);
                PureOutputMain(boolMat, idx, lenGrps, boolVec,
                               z, nCols, nRows, nThreads, Parallel);
                return res;
            } case RAWSXP : {
                cpp11::sexp res = Rf_allocMatrix(RAWSXP, nRows, nCols);
                Rbyte* rawMat = RAW(res);
                GetPureOutput(rawMat, idx, lenGrps, rawVec, z, nCols, nRows);
                return res;
            } case CPLXSXP : {
                cpp11::sexp res = Rf_allocMatrix(CPLXSXP, nRows, nCols);
                Rcomplex* cmplxMat = COMPLEX(res);
                GetPureOutput(cmplxMat, idx, lenGrps,
                              cmplxVec, z, nCols, nRows);
                return res;
            } case REALSXP : {
                cpp11::sexp res = Rf_allocMatrix(REALSXP, nRows, nCols);
                double* dblMat = REAL(res);
                PureOutputMain(dblMat, idx, lenGrps, dblVec,
                               z, nCols, nRows, nThreads, Parallel);
                return res;
            } case STRSXP : {
                cpp11::writable::strings_matrix<> charMat(nRows, nCols);
                GetCharOutput(charMat, idx, lenGrps,
                              charVec, z, nCols, nRows);
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

    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                 VecType::Integer, "maxThreads");

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

    // Transform lenGrps to be used in nextProduct
    for (auto &v_i: lenGrps) {
        v_i = nCols * (v_i - 1);
    }

    cpp11::sexp res = GlueProdCart(
        mat_idx, typeCheck, IsFactor, RList, intVec, dblVec,
        boolVec, cmplxVec, rawVec, charVec, lenGrps, startZ,
        nRows, nCols, IsDF, nThreads, Parallel
    );

    return res;
}
