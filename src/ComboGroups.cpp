#include "cpp11/strings.hpp"
#include "cpp11/list.hpp"

#include "ComboGroupsUtils.h"
#include "CppConvert.h"
#include "SetUpUtils.h"
#include "RMatrix.h"
#include <numeric>
#include <thread>

void FinalTouch(SEXP res, bool IsArray, int grpSize, int r, int n,
                int nRows, bool IsNamed, const std::vector<double> &mySample,
                const std::vector<mpz_class> &myBigSamp,
                bool IsGmp, bool IsSample) {

    std::vector<std::string> myColNames(r, "Grp");

    for (int j = 0; j < r; ++j) {
        myColNames[j] += std::to_string(j + 1);
    }

    if (IsArray) {
        cpp11::integers dim({nRows, grpSize, r});
        Rf_setAttrib(res, R_DimSymbol, dim);
        cpp11::writable::strings myNames(r);

        for (int i = 0; i < r; ++i) {
            myNames[i] = myColNames[i].c_str();
        }

        SetSampleNames(res, IsGmp, nRows, mySample,
                       myBigSamp, IsNamed, myNames, 2);

        if (!IsNamed) {
            cpp11::writable::list dimNames(3);
            dimNames[2] = myNames;
            Rf_setAttrib(res, R_DimNamesSymbol, dimNames);
        }
    } else {
        cpp11::writable::strings myNames(n);

        for (int i = 0, k = 0; i < r; ++i) {
            for (int j = 0; j < grpSize; ++j, ++k) {
                myNames[k] = myColNames[i].c_str();
            }
        }

        SetSampleNames(res, IsGmp, nRows, mySample,
                       myBigSamp, IsNamed, myNames, 1);

        if (!IsNamed) {
            cpp11::writable::list dimNames(2);
            dimNames[1] = myNames;
            Rf_setAttrib(res, R_DimNamesSymbol, dimNames);
        }
    }
}

void SampleResults(SEXP GroupsMat, SEXP v,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   const mpz_class &computedRowMpz, double computedRows,
                   int sampSize, int n, int r, int grpSize, bool IsGmp) {

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
                   const std::vector<mpz_class> &myBigSamp,
                   const mpz_class &computedRowMpz, double computedRows,
                   std::size_t sampSize, std::size_t n, int r,
                   int grpSize, bool IsGmp) {

    if (IsGmp) {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthComboGroupGmp(
                n, grpSize, r, myBigSamp[i], computedRowMpz
            );

            for (std::size_t j = 0; j < n; ++j) {
                GroupsMat[i + sampSize * j] = v[z[j]];
            }
        }
    } else {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> z = nthComboGroup(
                n, grpSize, r, mySample[i], computedRows
            );

            for (std::size_t j = 0; j < n; ++j) {
                GroupsMat[i + sampSize * j] = v[z[j]];
            }
        }
    }
}

template <typename T>
void SampleResults(RcppParallel::RMatrix<T> GroupsMat,
                   const std::vector<T> &v,
                   const std::vector<double> &mySample,
                   const std::vector<mpz_class> &myBigSamp,
                   const mpz_class &computedRowMpz, double computedRows,
                   std::size_t n, int r, int grpSize, bool IsGmp,
                   std::size_t strtIdx, std::size_t endIdx) {

    if (IsGmp) {
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthComboGroupGmp(
                n, grpSize, r, myBigSamp[i], computedRowMpz
            );

            for (std::size_t j = 0; j < n; ++j) {
                GroupsMat(i, j) = v[z[j]];
            }
        }
    } else {
        for (std::size_t i = strtIdx; i < endIdx; ++i) {
            const std::vector<int> z = nthComboGroup(
                n, grpSize, r, mySample[i], computedRows
            );

            for (std::size_t j = 0; j < n; ++j) {
                GroupsMat(i, j) = v[z[j]];
            }
        }
    }
}

void GroupWorker(SEXP GroupsMat, SEXP v, std::vector<int> &z,
                 std::size_t nRows, std::size_t n, int r, int grpSize) {

    const int idx1 = (r - 1) * grpSize - 1;
    const int idx2 = Rf_length(v) - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const std::size_t lastRow = nRows - 1;

    for (std::size_t i = 0; i < lastRow; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            SET_STRING_ELT(GroupsMat, i + j * nRows, STRING_ELT(v, z[j]));
        }

        nextComboGroup(z, r, grpSize, idx1, idx2, last1);
    }

    // Get last combo group
    for (std::size_t j = 0; j < n; ++j) {
        SET_STRING_ELT(GroupsMat, lastRow + j * nRows, STRING_ELT(v, z[j]));
    }
}

template <typename T>
void GroupWorker(T* GroupsMat, const std::vector<T> &v,
                 std::vector<int> &z, std::size_t nRows,
                 std::size_t n, int r, int grpSize) {

    const int idx1 = (r - 1) * grpSize - 1;
    const int idx2 = v.size() - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const std::size_t lastRow = nRows - 1;

    for (std::size_t i = 0; i < lastRow; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            GroupsMat[i + j * nRows] = v[z[j]];
        }

        nextComboGroup(z, r, grpSize, idx1, idx2, last1);
    }

    // Get last combo group
    for (std::size_t j = 0; j < n; ++j) {
        GroupsMat[lastRow + j * nRows] = v[z[j]];
    }
}

template <typename T>
void GroupWorker(RcppParallel::RMatrix<T> &GroupsMat,
                 const std::vector<T> &v, std::vector<int> &z, std::size_t n,
                 int r, int grpSize, std::size_t strtIdx, std::size_t endIdx) {

    const int idx1 = (r - 1) * grpSize - 1;
    const int idx2 = v.size() - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const std::size_t lastRow = endIdx - 1;

    for (std::size_t i = strtIdx; i < lastRow; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            GroupsMat(i, j) = v[z[j]];
        }

        nextComboGroup(z, r, grpSize, idx1, idx2, last1);
    }

    // Get last combo group
    for (std::size_t j = 0; j < n; ++j) {
        GroupsMat(lastRow, j) = v[z[j]];
    }
}

template <typename T>
void SerialGlue(T* GroupsMat, SEXP res, const std::vector<T> &v,
                const std::vector<double> &mySamp,
                const std::vector<mpz_class> &myBigSamp, std::vector<int> z,
                const mpz_class &computedRowMpz, double computedRows, int n,
                int r, int grpSize, int nRows, bool IsArray,
                bool IsGmp, bool IsSample, bool IsNamed) {

    if (IsSample) {
        SampleResults(GroupsMat, v, mySamp, myBigSamp, computedRowMpz,
                      computedRows, nRows, n, r, grpSize, IsGmp);
    } else {
        GroupWorker(GroupsMat, v, z, nRows, n, r, grpSize);
    }

    FinalTouch(res, IsArray, grpSize, r, n, nRows,
               IsNamed, mySamp, myBigSamp, IsGmp, IsSample);
}

void CharacterGlue(SEXP res, SEXP v,
                   const std::vector<double> &mySamp,
                   const std::vector<mpz_class> &myBigSamp, std::vector<int> z,
                   const mpz_class &computedRowMpz, double computedRows, int n,
                   int r, int grpSize, int nRows, bool IsArray,
                   bool IsGmp, bool IsSample, bool IsNamed) {

    if (IsSample) {
        SampleResults(res, v, mySamp, myBigSamp, computedRowMpz,
                      computedRows, nRows, n, r, grpSize, IsGmp);
    } else {
        GroupWorker(res, v, z, nRows, n, r, grpSize);
    }

    FinalTouch(res, IsArray, grpSize, r, n, nRows,
               IsNamed, mySamp, myBigSamp, IsGmp, IsSample);
}

template <typename T>
void ParallelGlue(RcppParallel::RMatrix<T> &GroupsMat,
                  const std::vector<T> &v, const std::vector<double> &mySamp,
                  const std::vector<mpz_class> &myBigSamp, std::vector<int> z,
                  const mpz_class &computedRowMpz, double computedRows, int n,
                  int r, int grpSize, int strtIdx, int endIdx,
                  bool IsGmp, bool IsSample) {

    if (IsSample) {
        SampleResults(GroupsMat, v, mySamp, myBigSamp, computedRowMpz,
                      computedRows, n, r, grpSize, IsGmp, strtIdx, endIdx);
    } else {
        GroupWorker(GroupsMat, v, z, n, r, grpSize, strtIdx, endIdx);
    }
}

void GetStartGrp(std::vector<int> &z,
                 const mpz_class &computedRowMpz, mpz_class &lowerMpz,
                 double computedRows, double &lower, int n, int grpSize,
                 int r, int stepSize, bool IsGmp) {

    if (IsGmp) {
        lowerMpz += stepSize;
        z = nthComboGroupGmp(n, grpSize, r, lowerMpz, computedRowMpz);
    } else {
        lower += stepSize;
        z = nthComboGroup(n, grpSize, r, lower, computedRows);
    }
}

template <typename T>
void GroupsMain(T* GroupsMat, SEXP res, const std::vector<T> &v,
                std::vector<int> z, const std::vector<double> &mySample,
                const std::vector<mpz_class> &myBigSamp,
                const mpz_class &computedRowMpz,
                double computedRows, mpz_class lowerMpz, double lower,
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
                std::cref(ParallelGlue<T>), std::ref(parMat),
                std::cref(v), std::cref(mySample), std::cref(myBigSamp),
                z, std::cref(computedRowMpz), computedRows, n, r,
                grpSize, step, nextStep, IsGmp, IsSample
            );

            GetStartGrp(z, computedRowMpz, lowerMpz, computedRows,
                        lower, n, grpSize, r, stepSize, IsGmp);
        }

        threads.emplace_back(std::cref(ParallelGlue<T>), std::ref(parMat),
                             std::cref(v), std::cref(mySample),
                             std::cref(myBigSamp), z, std::cref(computedRowMpz),
                             computedRows, n, r, grpSize, step, nRows,
                             IsGmp, IsSample);

        for (auto& thr: threads) {
            thr.join();
        }

        FinalTouch(res, IsArray, grpSize, r, n, nRows,
                   IsNamed, mySample, myBigSamp, IsGmp, IsSample);
    } else {
        SerialGlue(GroupsMat, res, v, mySample, myBigSamp, z, computedRowMpz,
                   computedRows, n, r, grpSize, nRows, IsArray, IsGmp,
                   IsSample, IsNamed);
    }
}

[[cpp11::register]]
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
    CppConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

    CppConvert::convertPrimitive(RNumGroups, numGroups,
                                   VecType::Integer, "numGroups");

    bool IsSample = CppConvert::convertFlag(RIsSample, "IsSample");
    bool Parallel = CppConvert::convertFlag(Rparallel, "Parallel");
    bool IsNamed = (IsSample) ? CppConvert::convertFlag(RNamed,
                    "namedSample") : false;

    std::vector<int> vInt;
    std::vector<double> vNum;

    SetType(myType, Rv);
    SetBasic(Rv, vNum, vInt, n, myType);

    if (myType == VecType::Integer) {
        vInt.assign(vNum.cbegin(), vNum.cend());
    }

    if (n % numGroups != 0) {
        cpp11::stop("The length of v (if v is a vector) or v (if v"
                 " is a scalar) must be divisible by numGroups");
    }

    const int grpSize = n / numGroups;
    const double computedRows = numGroupCombs(n, numGroups, grpSize);
    bool IsGmp = (computedRows > SampleLimit);
    mpz_class computedRowMpz;

    if (IsGmp) {
        computedRowMpz = 1;
        numGroupCombsGmp(computedRowMpz, n, numGroups, grpSize);
    } else {
        computedRowMpz = computedRows;
    }

    double lower = 0;
    double upper = 0;

    bool bLower = false;
    bool bUpper = false;

    mpz_class lowerMpz;
    mpz_class upperMpz;

    if (!IsSample) {
        SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
                  lowerMpz, upperMpz, computedRowMpz, computedRows);
    }

    std::vector<int> startZ;
    double dblLower = lower;
    if (!IsGmp) lowerMpz = dblLower;

    if (bLower && cmp(lowerMpz, 0) > 0) {
        startZ = IsGmp ? nthComboGroupGmp(n, grpSize, numGroups,
                                          lowerMpz, computedRowMpz) :
        nthComboGroup(n, grpSize, numGroups, lower, computedRows);
    } else {
        startZ.resize(n);
        std::iota(startZ.begin(), startZ.end(), 0);
    }

    int nRows = 0;

    if (!IsSample) {
        double userNumRows = 0;
        SetNumResults(IsGmp, bLower, bUpper, true, upperMpz,
                      lowerMpz, lower, upper, computedRows,
                      computedRowMpz, nRows, userNumRows);
    }

    const std::string retType(CHAR(STRING_ELT(RRetType, 0)));

    if (retType != "3Darray" && retType != "matrix") {
        cpp11::stop("retType must be '3Darray' or 'matrix'");
    }

    int sampSize;
    std::vector<double> mySample;
    const bool IsArray = (retType == "3Darray");

    if (IsSample) {
        SetRandomSample(RindexVec, RNumSamp, sampSize, IsGmp,
                        computedRows, mySample, baseSample, myEnv);
    }

    const int bigSampSize = (IsSample && IsGmp) ? sampSize : 1;
    std::vector<mpz_class> myVec(bigSampSize);

    if (IsSample) {
        SetRandomSampleMpz(RindexVec, RmySeed, sampSize,
                           IsGmp, computedRowMpz, myVec);
    }

    const int numResults = (IsSample) ? sampSize : nRows;
    const int limit = (IsSample) ? 2 : 20000;
    SetThreads(Parallel, maxThreads, numResults,
               myType, nThreads, RNumThreads, limit);

    switch (myType) {
        case VecType::Character: {
            cpp11::sexp charVec = Rf_duplicate(Rv);
            cpp11::sexp res = Rf_allocMatrix(STRSXP, numResults, n);

            CharacterGlue(res, charVec, mySample, myVec, startZ,
                          computedRowMpz, computedRows, n, numGroups, grpSize,
                          numResults, IsArray, IsGmp, IsSample, IsNamed);

            return res;
        } case VecType::Complex: {
            std::vector<Rcomplex> stlCmplxVec(n);
            Rcomplex* vecCmplx = COMPLEX(Rv);

            for (int i = 0; i < n; ++i) {
                stlCmplxVec[i] = vecCmplx[i];
            }

            cpp11::sexp res = Rf_allocMatrix(CPLXSXP, numResults, n);
            Rcomplex* matCmplx = COMPLEX(res);

            SerialGlue(matCmplx, res, stlCmplxVec, mySample, myVec,
                       startZ, computedRowMpz, computedRows, n, numGroups,
                       grpSize, numResults, IsArray, IsGmp, IsSample,
                       IsNamed);

            return res;
        } case VecType::Raw : {
            std::vector<Rbyte> stlRawVec(n);
            Rbyte* vecRaw = RAW(Rv);

            for (int i = 0; i < n; ++i) {
                stlRawVec[i] = vecRaw[i];
            }

            cpp11::sexp res = Rf_allocMatrix(RAWSXP, numResults, n);
            Rbyte* matRaw = RAW(res);

            SerialGlue(matRaw, res, stlRawVec, mySample, myVec, startZ,
                       computedRowMpz, computedRows, n, numGroups, grpSize,
                       numResults, IsArray, IsGmp, IsSample, IsNamed);

            return res;
        } case VecType::Logical : {
            vInt.assign(n, 0);
            int* vecBool = LOGICAL(Rv);

            for (int i = 0; i < n; ++i) {
                vInt[i] = vecBool[i];
            }

            cpp11::sexp res = Rf_allocMatrix(LGLSXP, numResults, n);
            int* matBool = LOGICAL(res);

            SerialGlue(matBool, res, vInt, mySample, myVec, startZ,
                       computedRowMpz, computedRows, n, numGroups, grpSize,
                       numResults, IsArray, IsGmp, IsSample, IsNamed);

            return res;
        } case VecType::Integer : {
            cpp11::sexp res = Rf_allocMatrix(INTSXP, numResults, n);
            int* matInt = INTEGER(res);

            GroupsMain(matInt, res, vInt, startZ, mySample, myVec,
                       computedRowMpz, computedRows, lowerMpz, lower, n,
                       numGroups, grpSize, numResults, nThreads, IsArray,
                       IsNamed, Parallel, IsGmp, IsSample);

            if (Rf_isFactor(Rv)) {
                SetFactorClass(res, Rv);
            }

            return res;
        } default : {
            cpp11::sexp res = Rf_allocMatrix(REALSXP, numResults, n);
            double* matNum = REAL(res);

            GroupsMain(matNum, res, vNum, startZ, mySample, myVec,
                       computedRowMpz, computedRows, lowerMpz, lower, n,
                       numGroups, grpSize, numResults, nThreads, IsArray,
                       IsNamed, Parallel, IsGmp, IsSample);

            return res;
        }
    }
}
