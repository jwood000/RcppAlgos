#include "Cartesian/GetProduct.h"

void SampleResults(
    cpp11::writable::strings_matrix<> &mat, const cpp11::strings &charVec,
    const std::vector<int> &idx, const std::vector<int> &lenNxtPr,
    const std::vector<double> &mySample,
    const std::vector<mpz_class> &myBigSamp, int sampSize, int n, bool IsGmp
) {

    if (IsGmp) {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthProductGmp(myBigSamp[i], lenNxtPr);

            for (int j = 0; j < n; ++j) {
                SET_STRING_ELT(mat, i + sampSize * j,
                               STRING_ELT(charVec, idx[j + z[j]]));
            }
        }
    } else {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthProduct(mySample[i], lenNxtPr);

            for (int j = 0; j < n; ++j) {
                SET_STRING_ELT(mat, i + sampSize * j,
                               STRING_ELT(charVec, idx[j + z[j]]));
            }
        }
    }
}

template <typename T>
void SampleResults(
    T* ProdMat, const std::vector<T> &v, const std::vector<int> &idx,
    const std::vector<int> &lenNxtPr, const std::vector<double> &mySample,
    const std::vector<mpz_class> &myBigSamp, std::size_t sampSize,
    std::size_t nCols, bool IsGmp
) {

    if (IsGmp) {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthProductGmp(myBigSamp[i], lenNxtPr);

            for (std::size_t j = 0; j < nCols; ++j) {
                ProdMat[i + sampSize * j] = v[idx[j + z[j]]];
            }
        }
    } else {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthProduct(mySample[i], lenNxtPr);

            for (std::size_t j = 0; j < nCols; ++j) {
                ProdMat[i + sampSize * j] = v[idx[j + z[j]]];
            }
        }
    }
}

template <typename T>
void SampleResults(
    RcppParallel::RMatrix<T> ProdMat, const std::vector<T> &v,
    const std::vector<int> &idx, const std::vector<int> &lenNxtPr,
    const std::vector<double> &mySample,
    const std::vector<mpz_class> &myBigSamp, std::size_t nCols,
    std::size_t strtIdx, std::size_t endIdx, bool IsGmp
) {

    if (IsGmp) {
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthProductGmp(myBigSamp[i], lenNxtPr);

            for (std::size_t j = 0; j < nCols; ++j) {
                ProdMat(i, j) = v[idx[j + z[j]]];
            }
        }
    } else {
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthProduct(mySample[i], lenNxtPr);

            for (std::size_t j = 0; j < nCols; ++j) {
                ProdMat(i, j) = v[idx[j + z[j]]];
            }
        }
    }
}

template <typename T>
void GetPureOutput(T* mat, const std::vector<int> &idx,
                   const std::vector<int> &lenGrps,
                   const std::vector<T> &v,
                   std::vector<int> &z, int nCols, int nRows) {

    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            mat[i + j * nRows] = v[idx[j + z[j]]];
        }

        nextProduct(lenGrps, z, nCols);
    }
}

template <typename T>
void ParallelProduct(
    RcppParallel::RMatrix<T> &mat, const std::vector<int> &idx,
    const std::vector<int> &lenGrps, const std::vector<T> &v,
    std::vector<int> z, int nCols, int strt, int nRows
) {

    for (int i = strt; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            mat(i, j) = v[idx[j + z[j]]];
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
void SerialGlue(
    T* ProdMat, const std::vector<int> &idx,
    const std::vector<int> &lenGrps, const std::vector<T> &v,
    const std::vector<int> &lenNxtPr, const std::vector<double> &mySamp,
    const std::vector<mpz_class> &myBigSamp,
    std::vector<int> z, int nCols, int nRows,
    bool IsSample, bool IsGmp
) {

    if (IsSample) {
        SampleResults(ProdMat, v, idx, lenNxtPr, mySamp,
                      myBigSamp, nRows, nCols, IsGmp);
    } else {
        GetPureOutput(ProdMat, idx, lenGrps, v, z, nCols, nRows);
    }
}

void CharacterGlue(
    cpp11::writable::strings_matrix<> &mat, const cpp11::strings &charVec,
    const std::vector<int> &idx, const std::vector<int> &lenGrps,
    const std::vector<int> &lenNxtPr, const std::vector<double> &mySamp,
    const std::vector<mpz_class> &myBigSamp,
    std::vector<int> z, int nCols, int nRows, bool IsSample, bool IsGmp
) {

    if (IsSample) {
        SampleResults(mat, charVec, idx, lenNxtPr, mySamp,
                      myBigSamp, nRows, nCols, IsGmp);
    } else {
        GetCharOutput(mat, idx, lenGrps, charVec, z, nCols, nRows);
    }
}

template <typename T>
void ParallelGlue(
    RcppParallel::RMatrix<T> &ProdMat, const std::vector<int> &idx,
    const std::vector<int> &lenGrps, const std::vector<T> &v,
    const std::vector<int> &lenNxtPr, const std::vector<double> &mySamp,
    const std::vector<mpz_class> &myBigSamp, std::vector<int> z,
    int nCols, int strt, int nRows, bool IsSample, bool IsGmp
) {

    if (IsSample) {
        SampleResults(ProdMat, v, idx, lenNxtPr, mySamp,
                      myBigSamp, nCols, strt, nRows, IsGmp);
    } else {
        ParallelProduct(
            ProdMat, idx, lenGrps, v, z, nCols, strt, nRows
        );
    }
}

template <typename T>
void PureOutputMain(
    T* mat, const std::vector<int> &idx,
    const std::vector<int> &lenGrps, const std::vector<T> &v,
    const std::vector<int> &lenNxtPr, const std::vector<double> &mySamp,
    const std::vector<mpz_class> &myBigSamp,
    std::vector<int> z, int nCols, int nRows, int nThreads, bool Parallel,
    mpz_class lowerMpz, double lower, bool IsSample, bool IsGmp
) {

    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(mat, nRows, nCols);
        std::vector<std::thread> threads;

        int step = 0;
        int stepSize = nRows / nThreads;
        int nextStep = stepSize;

        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize,
             nextStep += stepSize) {

            threads.emplace_back(
                std::cref(ParallelGlue<T>), std::ref(parMat),
                std::cref(idx), std::cref(lenGrps), std::cref(v),
                std::cref(lenNxtPr), std::cref(mySamp), std::cref(myBigSamp),
                z, nCols, step, nextStep, IsSample, IsGmp
            );

            GetStartProd(lenNxtPr, z, lowerMpz, lower, stepSize, IsGmp);
        }

        threads.emplace_back(
            std::cref(ParallelGlue<T>), std::ref(parMat),
            std::cref(idx), std::cref(lenGrps), std::cref(v),
            std::cref(lenNxtPr), std::cref(mySamp), std::cref(myBigSamp),
            z, nCols, step, nRows, IsSample, IsGmp
        );

        for (auto& thr: threads) {
            thr.join();
        }
    } else {
        SerialGlue(mat, idx, lenGrps, v, lenNxtPr, mySamp,
                   myBigSamp, z, nCols, nRows, IsSample, IsGmp);
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

SEXP GetProduct(
    const std::vector<int> &idx, const std::vector<int> &typeCheck,
    const std::vector<int> &IsFactor, const cpp11::list &RList,
    const std::vector<int> &intVec, const std::vector<double> &dblVec,
    const std::vector<int> &boolVec, const std::vector<Rcomplex> &cmplxVec,
    const std::vector<Rbyte> &rawVec, const cpp11::strings &charVec,
    const std::vector<int> &lenGrps, std::vector<int> &z,
    const std::vector<double> &mySamp, const std::vector<mpz_class> &myBigSamp,
    double lower, mpz_class &lowerMpz, int nRows, int nCols, bool IsDF,
    int nThreads, bool Parallel, bool IsGmp, bool IsSample
) {

    std::vector<int> lenNxtPr(lenGrps);

    for (auto &v_i: lenNxtPr) {
        v_i = (v_i / nCols) + 1;
    }

    if (IsDF) {
        cpp11::writable::list DataFrame(nCols);
        std::vector<int> all_idx(nRows * nCols);

        if (IsSample && IsGmp) {
            for (int i = 0; i < nRows; ++i) {
                const std::vector<int> z =
                    nthProductGmp(myBigSamp[i], lenNxtPr);
                std::copy(z.begin(), z.end(), all_idx.begin() + i * nCols);
            }
        } else if (IsSample) {
            for (int i = 0; i < nRows; ++i) {
                const std::vector<int> z = nthProduct(mySamp[i], lenNxtPr);
                std::copy(z.begin(), z.end(), all_idx.begin() + i * nCols);
            }
        } else {
            for (int i = 0; i < nRows; ++i) {
                std::copy(z.begin(), z.end(), all_idx.begin() + i * nCols);
                nextProduct(lenGrps, z, nCols);
            }
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

                PureOutputMain(
                    intMat, idx, lenGrps, intVec, lenNxtPr, mySamp,
                    myBigSamp, z, nCols, nRows, nThreads, Parallel,
                    lowerMpz, lower, IsSample, IsGmp
                );

                if (typeCheck[tFac]) SetFactorClass(res, RList[0]);
                return res;
            } case LGLSXP : {
                cpp11::sexp res = Rf_allocMatrix(LGLSXP, nRows, nCols);
                int* boolMat = LOGICAL(res);

                PureOutputMain(
                    boolMat, idx, lenGrps, boolVec, lenNxtPr, mySamp,
                    myBigSamp, z, nCols, nRows, nThreads, Parallel,
                    lowerMpz, lower, IsSample, IsGmp
                );

                return res;
            } case RAWSXP : {
                cpp11::sexp res = Rf_allocMatrix(RAWSXP, nRows, nCols);
                Rbyte* rawMat = RAW(res);

                SerialGlue(rawMat, idx, lenGrps, rawVec, lenNxtPr, mySamp,
                           myBigSamp, z, nCols, nRows, IsSample, IsGmp);

                return res;
            } case CPLXSXP : {
                cpp11::sexp res = Rf_allocMatrix(CPLXSXP, nRows, nCols);
                Rcomplex* cmplxMat = COMPLEX(res);

                SerialGlue(cmplxMat, idx, lenGrps, cmplxVec, lenNxtPr, mySamp,
                           myBigSamp, z, nCols, nRows, IsSample, IsGmp);

                return res;
            } case REALSXP : {
                cpp11::sexp res = Rf_allocMatrix(REALSXP, nRows, nCols);
                double* dblMat = REAL(res);

                PureOutputMain(
                    dblMat, idx, lenGrps, dblVec, lenNxtPr, mySamp,
                    myBigSamp, z, nCols, nRows, nThreads, Parallel,
                    lowerMpz, lower, IsSample, IsGmp
                );

                return res;
            } case STRSXP : {
                cpp11::writable::strings_matrix<> charMat(nRows, nCols);

                CharacterGlue(charMat, charVec, idx, lenGrps, lenNxtPr, mySamp,
                              myBigSamp, z, nCols, nRows, IsSample, IsGmp);

                return charMat;
            } default : {
                cpp11::stop("Only atomic types are supported for v");
            }
        }
    }
}
