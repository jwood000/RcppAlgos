#include "ComboGroups/GetComboGroups.h"

void SampleResults(SEXP GroupsMat, SEXP v,
                   nthFuncDbl nthCmbGrp, nthFuncGmp nthCmbGrpGmp,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   int sampSize, int n, bool IsGmp) {

    if (IsGmp) {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthCmbGrpGmp(myBigSamp[i]);

            for (int j = 0; j < n; ++j) {
                SET_STRING_ELT(GroupsMat, i + sampSize * j,
                               STRING_ELT(v, z[j]));
            }
        }
    } else {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthCmbGrp(mySample[i]);

            for (int j = 0; j < n; ++j) {
                SET_STRING_ELT(GroupsMat, i + sampSize * j,
                               STRING_ELT(v, z[j]));
            }
        }
    }
}

template <typename T>
void SampleResults(T* GroupsMat, const std::vector<T> &v,
                   nthFuncDbl nthCmbGrp, nthFuncGmp nthCmbGrpGmp,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   std::size_t sampSize, std::size_t n, bool IsGmp) {

    if (IsGmp) {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthCmbGrpGmp(myBigSamp[i]);

            for (std::size_t j = 0; j < n; ++j) {
                GroupsMat[i + sampSize * j] = v[z[j]];
            }
        }
    } else {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthCmbGrp(mySample[i]);

            for (std::size_t j = 0; j < n; ++j) {
                GroupsMat[i + sampSize * j] = v[z[j]];
            }
        }
    }
}

template <typename T>
void SampleResults(RcppParallel::RMatrix<T> GroupsMat, const std::vector<T> &v,
                   nthFuncDbl nthCmbGrp, nthFuncGmp nthCmbGrpGmp,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp, std::size_t n,
                   std::size_t strtIdx, std::size_t endIdx, bool IsGmp) {

    if (IsGmp) {
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthCmbGrpGmp(myBigSamp[i]);

            for (std::size_t j = 0; j < n; ++j) {
                GroupsMat(i, j) = v[z[j]];
            }
        }
    } else {
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthCmbGrp(mySample[i]);

            for (std::size_t j = 0; j < n; ++j) {
                GroupsMat(i, j) = v[z[j]];
            }
        }
    }
}

void GroupWorker(SEXP GroupsMat, SEXP v, nextGrpFunc nextCmbGrp,
                 std::vector<int> &z, std::size_t nRows, std::size_t n) {

    const std::size_t lastRow = nRows - 1;

    for (std::size_t i = 0; i < lastRow; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            SET_STRING_ELT(GroupsMat, i + j * nRows, STRING_ELT(v, z[j]));
        }

        nextCmbGrp(z);
    }

    // Get last combo group
    for (std::size_t j = 0; j < n; ++j) {
        SET_STRING_ELT(GroupsMat, lastRow + j * nRows, STRING_ELT(v, z[j]));
    }
}

template <typename T>
void GroupWorker(T* GroupsMat, const std::vector<T> &v, nextGrpFunc nextCmbGrp,
                 std::vector<int> &z, std::size_t nRows, std::size_t n) {

    const std::size_t lastRow = nRows - 1;

    for (std::size_t i = 0; i < lastRow; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            GroupsMat[i + j * nRows] = v[z[j]];
        }

        nextCmbGrp(z);
    }

    // Get last combo group
    for (std::size_t j = 0; j < n; ++j) {
        GroupsMat[lastRow + j * nRows] = v[z[j]];
    }
}

template <typename T>
void GroupWorker(RcppParallel::RMatrix<T> &GroupsMat, const std::vector<T> &v,
                 nextGrpFunc nextCmbGrp, std::vector<int> &z, std::size_t n,
                 std::size_t strtIdx, std::size_t endIdx) {

    const std::size_t lastRow = endIdx - 1;

    for (std::size_t i = strtIdx; i < lastRow; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            GroupsMat(i, j) = v[z[j]];
        }

        nextCmbGrp(z);
    }

    // Get last combo group
    for (std::size_t j = 0; j < n; ++j) {
        GroupsMat(lastRow, j) = v[z[j]];
    }
}

template <typename T>
void SerialGlue(T* GroupsMat, SEXP res, const std::vector<T> &v,
                nextGrpFunc nextCmbGrp, nthFuncDbl nthCmbGrp,
                nthFuncGmp nthCmbGrpGmp, finalTouchFunc FinalTouch,
                const std::vector<double> &mySamp,
                const std::vector<mpz_class> &myBigSamp,
                std::vector<int> z, int n, int nRows, bool IsArray,
                bool IsSample, bool IsNamed, bool IsGmp) {

    if (IsSample) {
        SampleResults(GroupsMat, v, nthCmbGrp, nthCmbGrpGmp,
                      mySamp, myBigSamp, nRows, n, IsGmp);
    } else {
        GroupWorker(GroupsMat, v, nextCmbGrp, z, nRows, n);
    }

    FinalTouch(res, IsArray, nRows, IsNamed, mySamp, myBigSamp, IsSample);
}

void CharacterGlue(SEXP res, SEXP v, nextGrpFunc nextCmbGrp,
                   nthFuncDbl nthCmbGrp, nthFuncGmp nthCmbGrpGmp,
                   finalTouchFunc FinalTouch,
                   const std::vector<double> &mySamp,
                   const std::vector<mpz_class> &myBigSamp,
                   std::vector<int> z, int n, int nRows, bool IsArray,
                   bool IsSample, bool IsNamed, bool IsGmp) {

    if (IsSample) {
        SampleResults(res, v, nthCmbGrp, nthCmbGrpGmp,
                      mySamp, myBigSamp, nRows, n, IsGmp);
    } else {
        GroupWorker(res, v, nextCmbGrp, z, nRows, n);
    }

    FinalTouch(res, IsArray, nRows, IsNamed, mySamp, myBigSamp, IsSample);
}

template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T> &GroupsMat, const std::vector<T> &v,
                  nextGrpFunc nextCmbGrp, nthFuncDbl nthCmbGrp,
                  nthFuncGmp nthCmbGrpGmp, const std::vector<double> &mySamp,
                  const std::vector<mpz_class> &myBigSamp, std::vector<int> z,
                  int n, int strtIdx, int endIdx, bool IsSample, bool IsGmp) {

    if (IsSample) {
        SampleResults(GroupsMat, v, nthCmbGrp, nthCmbGrpGmp,
                      mySamp, myBigSamp, n, strtIdx, endIdx, IsGmp);
    } else {
        GroupWorker(GroupsMat, v, nextCmbGrp, z, n, strtIdx, endIdx);
    }
}

void GetStartGrp(nthFuncDbl nthCmbGrp, nthFuncGmp nthCmbGrpGmp,
                 std::vector<int> &z, mpz_class &lowerMpz,
                 double &lower, int stepSize, bool IsGmp) {

    if (IsGmp) {
        lowerMpz += stepSize;
        z = nthCmbGrpGmp(lowerMpz);
    } else {
        lower += stepSize;
        z = nthCmbGrp(lower);
    }
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
    std::vector<int> &z, int nCols, int nRows, int nThreads, bool Parallel,
    mpz_class lowerMpz, double lower
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

            if (IsGmp) {
                lowerMpz += nextStep;
                z = nthProductGmp(lowerMpz, lenNxtPr);
            } else {
                lower += nextStep;
                z = nthProduct(lower, lenNxtPr);
            }
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
                } default : {
                    cpp11::stop("Only atomic types are supported for v");
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
            } default : {
                cpp11::stop("Only atomic types are supported for v");
            }
        }
    }
}
