#include "ComboGroupsUtils.h"
#include "Cpp14MakeUnique.h"
#include "CleanConvert.h"
#include "ComboGroups.h"
#include "SetUpUtils.h"
#include "RMatrix.h"
#include <numeric>
#include <thread>

void FinalTouch(SEXP res, bool IsArray, int grpSize, int r, int n,
                int nRows, bool IsNamed, const std::vector<double> &mySample,
                mpz_t *const myBigSamp, bool IsGmp) {

    std::vector<std::string> myColNames(r, "Grp");

    for (int j = 0; j < r; ++j) {
        myColNames[j] += std::to_string(j + 1);
    }

    if (IsArray) {
        SEXP dim = PROTECT(Rf_allocVector(INTSXP, 3));
        INTEGER(dim)[0] = nRows;
        INTEGER(dim)[1] = grpSize;
        INTEGER(dim)[2] = r;
        Rf_setAttrib(res, R_DimSymbol, dim);
        UNPROTECT(1);

        SEXP myNames = PROTECT(Rf_allocVector(STRSXP, r));

        for (int i = 0; i < r; ++i) {
            SET_STRING_ELT(myNames, i, Rf_mkChar(myColNames[i].c_str()));
        }

        if (IsNamed) {
            SetSampleNames(res, IsGmp, nRows, mySample, myBigSamp, myNames, 2);
        } else {
            SEXP dimNames = PROTECT(Rf_allocVector(VECSXP, 3));
            SET_VECTOR_ELT(dimNames, 2, myNames);
            Rf_setAttrib(res, R_DimNamesSymbol, dimNames);
            UNPROTECT(2);
        }
    } else {
        SEXP myNames = PROTECT(Rf_allocVector(STRSXP, n));

        for (int i = 0, k = 0; i < r; ++i) {
            for (int j = 0; j < grpSize; ++j, ++k) {
                SET_STRING_ELT(myNames, k, Rf_mkChar(myColNames[i].c_str()));
            }
        }

        if (IsNamed) {
            SetSampleNames(res, IsGmp, nRows, mySample, myBigSamp, myNames, 1);
        } else {
            SEXP dimNames = PROTECT(Rf_allocVector(VECSXP, 2));
            SET_VECTOR_ELT(dimNames, 1, myNames);
            Rf_setAttrib(res, R_DimNamesSymbol, dimNames);
            UNPROTECT(2);
        }
    }
}

void SampleResults(SEXP GroupsMat, SEXP v,
                   const std::vector<double> &mySample,
                   mpz_t *const myBigSamp, mpz_t computedRowMpz,
                   double computedRows, int sampSize, int n,
                   int r, int grpSize, bool IsGmp) {

    if (IsGmp) {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthComboGroupGmp(
                n, grpSize, r, myBigSamp[i], computedRowMpz
            );

            for (int j = 0; j < n; ++j) {
                SET_STRING_ELT(GroupsMat, i + sampSize * j,
                               STRING_ELT(v, z[j]));
            }
        }
    } else {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthComboGroup(
                n, grpSize, r, mySample[i], computedRows
            );

            for (int j = 0; j < n; ++j) {
                SET_STRING_ELT(GroupsMat, i + sampSize * j,
                               STRING_ELT(v, z[j]));
            }
        }
    }
}

template <typename T>
void SampleResults(T* GroupsMat, const std::vector<T> &v,
                   const std::vector<double> &mySample,
                   mpz_t *const myBigSamp, mpz_t computedRowMpz,
                   double computedRows, int sampSize, int n,
                   int r, int grpSize, bool IsGmp) {

    if (IsGmp) {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthComboGroupGmp(
                n, grpSize, r, myBigSamp[i], computedRowMpz
            );

            for (int j = 0; j < n; ++j) {
                GroupsMat[i + sampSize * j] = v[z[j]];
            }
        }
    } else {
        for (int i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthComboGroup(
                n, grpSize, r, mySample[i], computedRows
            );

            for (int j = 0; j < n; ++j) {
                GroupsMat[i + sampSize * j] = v[z[j]];
            }
        }
    }
}

template <typename T>
void SampleResults(RcppParallel::RMatrix<T> GroupsMat,
                   const std::vector<T> &v,
                   const std::vector<double> &mySample,
                   mpz_t *const myBigSamp, mpz_t computedRowMpz,
                   double computedRows, int n, int r, int grpSize,
                   bool IsGmp, int strtIdx, int endIdx) {

    if (IsGmp) {
        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthComboGroupGmp(
                n, grpSize, r, myBigSamp[i], computedRowMpz
            );

            for (int j = 0; j < n; ++j) {
                GroupsMat(i, j) = v[z[j]];
            }
        }
    } else {
        for (int i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthComboGroup(
                n, grpSize, r, mySample[i], computedRows
            );

            for (int j = 0; j < n; ++j) {
                GroupsMat(i, j) = v[z[j]];
            }
        }
    }
}

void GroupWorker(SEXP GroupsMat, SEXP v, std::vector<int> &z,
                 int nRows, int n, int r, int grpSize) {

    const int idx1 = (r - 1) * grpSize - 1;
    const int idx2 = Rf_length(v) - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const int lastRow = nRows - 1;

    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < n; ++j) {
            SET_STRING_ELT(GroupsMat, i + j * nRows, STRING_ELT(v, z[j]));
        }

        nextComboGroup(z, r, grpSize, idx1, idx2, last1);
    }

    // Get last combo group
    for (int j = 0; j < n; ++j) {
        SET_STRING_ELT(GroupsMat, lastRow + j * nRows, STRING_ELT(v, z[j]));
    }
}

template <typename T>
void GroupWorker(T* GroupsMat, const std::vector<T> &v,
                 std::vector<int> &z, int nRows, int n,
                 int r, int grpSize) {

    const int idx1 = (r - 1) * grpSize - 1;
    const int idx2 = v.size() - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const int lastRow = nRows - 1;

    for (int i = 0; i < lastRow; ++i) {
        for (int j = 0; j < n; ++j) {
            GroupsMat[i + j * nRows] = v[z[j]];
        }

        nextComboGroup(z, r, grpSize, idx1, idx2, last1);
    }

    // Get last combo group
    for (int j = 0; j < n; ++j) {
        GroupsMat[lastRow + j * nRows] = v[z[j]];
    }
}

template <typename T>
void GroupWorker(RcppParallel::RMatrix<T> &GroupsMat,
                 const std::vector<T> &v, std::vector<int> &z, int n,
                 int r, int grpSize, int strtIdx, int endIdx) {

    const int idx1 = (r - 1) * grpSize - 1;
    const int idx2 = v.size() - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const int lastRow = endIdx - 1;

    for (int i = strtIdx; i < lastRow; ++i) {
        for (int j = 0; j < n; ++j) {
            GroupsMat(i, j) = v[z[j]];
        }

        nextComboGroup(z, r, grpSize, idx1, idx2, last1);
    }

    // Get last combo group
    for (int j = 0; j < n; ++j) {
        GroupsMat(lastRow, j) = v[z[j]];
    }
}

template <typename T>
void SerialGlue(T* GroupsMat, SEXP res, const std::vector<T> &v,
                const std::vector<double> &mySamp,
                mpz_t *const myBigSamp, std::vector<int> z,
                mpz_t computedRowMpz, double computedRows, int n,
                int r, int grpSize, int nRows, bool IsArray,
                bool IsGmp, bool IsSample, bool IsNamed) {

    if (IsSample) {
        SampleResults(GroupsMat, v, mySamp, myBigSamp, computedRowMpz,
                      computedRows, nRows, n, r, grpSize, IsGmp);
    } else {
        GroupWorker(GroupsMat, v, z, nRows, n, r, grpSize);
    }

    FinalTouch(res, IsArray, grpSize, r, n, nRows,
               IsNamed, mySamp, myBigSamp, IsGmp);
}

void CharacterGlue(SEXP res, SEXP v,
                   const std::vector<double> &mySamp,
                   mpz_t *const myBigSamp, std::vector<int> z,
                   mpz_t computedRowMpz, double computedRows, int n,
                   int r, int grpSize, int nRows, bool IsArray,
                   bool IsGmp, bool IsSample, bool IsNamed) {

    if (IsSample) {
        SampleResults(res, v, mySamp, myBigSamp, computedRowMpz,
                      computedRows, nRows, n, r, grpSize, IsGmp);
    } else {
        GroupWorker(res, v, z, nRows, n, r, grpSize);
    }

    FinalTouch(res, IsArray, grpSize, r, n, nRows,
               IsNamed, mySamp, myBigSamp, IsGmp);
}

template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T> &GroupsMat,
                  const std::vector<T> &v, const std::vector<double> &mySamp,
                  mpz_t *const myBigSamp, std::vector<int> z,
                  mpz_t computedRowMpz, double computedRows, int n,
                  int r, int grpSize, int strtIdx, int endIdx,
                  bool IsGmp, bool IsSample) {

    if (IsSample) {
        SampleResults(GroupsMat, v, mySamp, myBigSamp, computedRowMpz,
                      computedRows, n, r, grpSize, IsGmp, strtIdx, endIdx);
    } else {
        GroupWorker(GroupsMat, v, z, n, r, grpSize, strtIdx, endIdx);
    }
}

void GetStartGrp(std::vector<int> &z, mpz_t computedRowMpz, mpz_t lowerMpz,
                 double computedRows, double &lower, int n, int grpSize,
                 int r, int stepSize, bool IsGmp) {

    if (IsGmp) {
        mpz_add_ui(lowerMpz, lowerMpz, stepSize);
        z = nthComboGroupGmp(n, grpSize, r, lowerMpz, computedRowMpz);
    } else {
        lower += stepSize;
        z = nthComboGroup(n, grpSize, r, lower, computedRows);
    }
}

template <typename T>
void GroupsMain(T* GroupsMat, SEXP res, const std::vector<T> &v,
                std::vector<int> z, const std::vector<double> &mySample,
                mpz_t *const myBigSamp, mpz_t computedRowMpz,
                double computedRows, mpz_t lowerMpz, double lower,
                int n, int r, int grpSize, int nRows, int nThreads,
                bool IsArray, bool IsNamed, bool Parallel,
                bool IsGmp,  bool IsSample) {

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
                std::cref(mySample), myBigSamp, z, computedRowMpz,
                computedRows, n, r, grpSize, step, nextStep, IsGmp, IsSample
            );

            GetStartGrp(z, computedRowMpz, lowerMpz, computedRows,
                        lower, n, grpSize, r, stepSize, IsGmp);
        }

        threads.emplace_back(std::cref(ParallelGlue<T>), std::ref(parMat),
                             std::cref(v), std::cref(mySample), myBigSamp,
                             z, computedRowMpz, computedRows, n, r,
                             grpSize, step, nRows, IsGmp, IsSample);

        for (auto& thr: threads) {
            thr.join();
        }

        FinalTouch(res, IsArray, grpSize, r, n, nRows,
                   IsNamed, mySample, myBigSamp, IsGmp);
    } else {
        SerialGlue(GroupsMat, res, v, mySample, myBigSamp, z, computedRowMpz,
                   computedRows, n, r, grpSize, nRows, IsArray, IsGmp,
                   IsSample, IsNamed);
    }
}

SEXP ComboGroupsCpp(SEXP Rv, SEXP RNumGroups, SEXP RRetType, SEXP Rlow,
                    SEXP Rhigh, SEXP Rparallel, SEXP RNumThreads,
                    SEXP RmaxThreads, SEXP RIsSample, SEXP RindexVec,
                    SEXP RmySeed, SEXP RNumSamp, SEXP baseSample,
                    SEXP RNamed, SEXP myEnv) {

    int n;
    int numGroups;
    int nThreads = 1;
    int maxThreads = 1;

    VecType myType = VecType::Integer;
    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

    CleanConvert::convertPrimitive(RNumGroups, numGroups,
                                   VecType::Integer, "numGroups");

    bool IsSample = CleanConvert::convertFlag(RIsSample, "IsSample");
    bool Parallel = CleanConvert::convertFlag(Rparallel, "Parallel");
    bool IsNamed = (IsSample) ? CleanConvert::convertFlag(RNamed,
                    "namedSample") : false;

    std::vector<int> vInt;
    std::vector<double> vNum;

    SetType(myType, Rv);
    SetBasic(Rv, vNum, vInt, n, myType);

    if (myType == VecType::Integer) {
        vInt.assign(vNum.cbegin(), vNum.cend());
    }

    if (n % numGroups != 0) {
        Rf_error("The length of v (if v is a vector) or v (if v"
                 " is a scalar) must be divisible by numGroups");
    }

    const int grpSize = n / numGroups;
    const double computedRows = numGroupCombs(n, numGroups, grpSize);
    bool IsGmp = (computedRows > SampleLimit);

    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);

    if (IsGmp) {
        mpz_set_ui(computedRowMpz, 1);
        numGroupCombsGmp(computedRowMpz, n, numGroups, grpSize);
    } else {
        mpz_set_d(computedRowMpz, computedRows);
    }

    double lower = 0;
    double upper = 0;

    bool bLower = false;
    bool bUpper = false;

    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);

    mpz_init(lowerMpz[0]);
    mpz_init(upperMpz[0]);

    if (!IsSample) {
        SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
                  lowerMpz.get(), upperMpz.get(), computedRowMpz, computedRows);
    }

    std::vector<int> startZ;
    double dblLower = lower;
    if (!IsGmp) mpz_set_d(lowerMpz[0], dblLower);

    if (bLower && mpz_cmp_ui(lowerMpz[0], 0) > 0) {
        startZ = IsGmp ? nthComboGroupGmp(n, grpSize, numGroups,
                                          lowerMpz[0], computedRowMpz) :
        nthComboGroup(n, grpSize, numGroups, lower, computedRows);
    } else {
        startZ.resize(n);
        std::iota(startZ.begin(), startZ.end(), 0);
    }

    int nRows = 0;

    if (!IsSample) {
        double userNumRows = 0;
        SetNumResults(IsGmp, bLower, bUpper, true, upperMpz[0],
                      lowerMpz[0], lower, upper, computedRows,
                      computedRowMpz, nRows, userNumRows);
    }

    const std::string retType(CHAR(STRING_ELT(RRetType, 0)));

    if (retType != "3Darray" && retType != "matrix") {
        Rf_error("retType must be '3Darray' or 'matrix'");
    }

    int sampSize;
    std::vector<double> mySample;
    const bool IsArray = (retType == "3Darray");

    if (IsSample) {
        SetRandomSample(RindexVec, RNumSamp, sampSize, IsGmp,
                        computedRows, mySample, baseSample, myEnv);
    }

    const int bigSampSize = (IsSample && IsGmp) ? sampSize : 1;
    auto myVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);

    for (int i = 0; i < bigSampSize; ++i) {
        mpz_init(myVec[i]);
    }

    if (IsSample) {
        SetRandomSampleMpz(RindexVec, RmySeed, sampSize,
                           IsGmp, computedRowMpz, myVec.get());
    }

    const int numResults = (IsSample) ? sampSize : nRows;
    const int limit = (IsSample) ? 2 : 20000;
    SetThreads(Parallel, maxThreads, numResults,
               myType, nThreads, RNumThreads, limit);

    switch (myType) {
        case VecType::Character: {
            SEXP charVec = PROTECT(Rf_duplicate(Rv));
            SEXP res = PROTECT(Rf_allocMatrix(STRSXP, numResults, n));

            CharacterGlue(res, charVec, mySample, myVec.get(), startZ,
                          computedRowMpz, computedRows, n, numGroups, grpSize,
                          numResults, IsArray, IsGmp, IsSample, IsNamed);

            UNPROTECT(2);
            return res;
        } case VecType::Complex: {
            std::vector<Rcomplex> stlCmplxVec(n);
            Rcomplex* vecCmplx = COMPLEX(Rv);

            for (int i = 0; i < n; ++i) {
                stlCmplxVec[i] = vecCmplx[i];
            }

            SEXP res = PROTECT(Rf_allocMatrix(CPLXSXP, numResults, n));
            Rcomplex* matCmplx = COMPLEX(res);

            SerialGlue(matCmplx, res, stlCmplxVec, mySample, myVec.get(),
                       startZ, computedRowMpz, computedRows, n, numGroups,
                       grpSize, numResults, IsArray, IsGmp, IsSample,
                       IsNamed);

            UNPROTECT(1);
            return res;
        } case VecType::Raw : {
            std::vector<Rbyte> stlRawVec(n);
            Rbyte* vecRaw = RAW(Rv);

            for (int i = 0; i < n; ++i) {
                stlRawVec[i] = vecRaw[i];
            }

            SEXP res = PROTECT(Rf_allocMatrix(RAWSXP, numResults, n));
            Rbyte* matRaw = RAW(res);

            SerialGlue(matRaw, res, stlRawVec, mySample, myVec.get(), startZ,
                       computedRowMpz, computedRows, n, numGroups, grpSize,
                       numResults, IsArray, IsGmp, IsSample, IsNamed);

            UNPROTECT(1);
            return res;
        } case VecType::Logical : {
            vInt.assign(n, 0);
            int* vecBool = LOGICAL(Rv);

            for (int i = 0; i < n; ++i) {
                vInt[i] = vecBool[i];
            }

            SEXP res = PROTECT(Rf_allocMatrix(LGLSXP, numResults, n));
            int* matBool = LOGICAL(res);

            SerialGlue(matBool, res, vInt, mySample, myVec.get(), startZ,
                       computedRowMpz, computedRows, n, numGroups, grpSize,
                       numResults, IsArray, IsGmp, IsSample, IsNamed);

            UNPROTECT(1);
            return res;
        } case VecType::Integer : {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, numResults, n));
            int* matInt = INTEGER(res);

            GroupsMain(matInt, res, vInt, startZ, mySample, myVec.get(),
                       computedRowMpz, computedRows, lowerMpz[0], lower, n,
                       numGroups, grpSize, numResults, nThreads, IsArray,
                       IsNamed, Parallel, IsGmp, IsSample);

            if (Rf_isFactor(Rv)) {
                SetFactorClass(res, Rv);
            }

            UNPROTECT(1);
            return res;
        } default : {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, numResults, n));
            double* matNum = REAL(res);

            GroupsMain(matNum, res, vNum, startZ, mySample, myVec.get(),
                       computedRowMpz, computedRows, lowerMpz[0], lower, n,
                       numGroups, grpSize, numResults, nThreads, IsArray,
                       IsNamed, Parallel, IsGmp, IsSample);

            UNPROTECT(1);
            return res;
        }
    }
}
