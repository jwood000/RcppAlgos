#include "ComboGroups/GetComboGroups.h"

CmbGrpClsFuncs GetClassFuncs(
    std::unique_ptr<ComboGroupsTemplate> const &CmbGrp
) {

    const nextGrpFunc nextCmbGrp = std::bind(
        &ComboGroupsTemplate::nextComboGroup,
        std::cref(CmbGrp), std::placeholders::_1
    );

    const nthFuncDbl nthCmbGrp = std::bind(
        &ComboGroupsTemplate::nthComboGroup,
        std::cref(CmbGrp), std::placeholders::_1
    );

    const nthFuncGmp nthCmbGrpGmp = std::bind(
        &ComboGroupsTemplate::nthComboGroupGmp,
        std::cref(CmbGrp), std::placeholders::_1
    );

    const finalTouchFunc FinalTouch = std::bind(
        &ComboGroupsTemplate::FinalTouch, std::cref(CmbGrp),
        std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
        std::placeholders::_4, std::placeholders::_5, std::placeholders::_6,
        std::placeholders::_7
    );

    CmbGrpClsFuncs myFuncs = {
        nthCmbGrp,
        nthCmbGrpGmp,
        nextCmbGrp,
        FinalTouch
    };

    return myFuncs;
}

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

SEXP GetComboGroups(
    SEXP Rv, nextGrpFunc nextCmbGrp, nthFuncDbl nthCmbGrp,
    nthFuncGmp nthCmbGrpGmp, finalTouchFunc FinalTouch,
    const std::vector<double> &vNum, const std::vector<int> &vInt,
    const std::vector<int> &startZ, const VecType &myType,
    const std::vector<double> &mySample,
    const std::vector<mpz_class> &myBigSamp,
    mpz_class lowerMpz, double lower, int n, int numResults, int nThreads,
    bool IsArray, bool IsNamed, bool Parallel, bool IsSample, bool IsGmp
) {

    switch (myType) {
        case VecType::Character: {
            cpp11::sexp charVec = Rf_duplicate(Rv);
            cpp11::sexp res = Rf_allocMatrix(STRSXP, numResults, n);

            CharacterGlue(res, charVec, nextCmbGrp, nthCmbGrp, nthCmbGrpGmp,
                          FinalTouch, mySample, myBigSamp, startZ, n,
                          numResults, IsArray, IsSample, IsNamed, IsGmp);

            return res;
        } case VecType::Complex: {
            std::vector<Rcomplex> stlCpxVec = CppConvert::GetVec<Rcomplex>(Rv);
            cpp11::sexp res = Rf_allocMatrix(CPLXSXP, numResults, n);
            Rcomplex* matCmplx = COMPLEX(res);

            SerialGlue(matCmplx, res, stlCpxVec, nextCmbGrp, nthCmbGrp,
                       nthCmbGrpGmp, FinalTouch, mySample, myBigSamp, startZ,
                       n, numResults, IsArray, IsSample, IsNamed, IsGmp);

            return res;
        } case VecType::Raw : {
            std::vector<Rbyte> stlRawVec = CppConvert::GetVec<Rbyte>(Rv);
            cpp11::sexp res = Rf_allocMatrix(RAWSXP, numResults, n);
            Rbyte* matRaw = RAW(res);

            SerialGlue(matRaw, res, stlRawVec, nextCmbGrp, nthCmbGrp,
                       nthCmbGrpGmp, FinalTouch, mySample, myBigSamp, startZ,
                       n, numResults, IsArray, IsSample, IsNamed, IsGmp);

            return res;
        } case VecType::Logical : {
            std::vector<int> vIntBool(n);
            int* vecBool = LOGICAL(Rv);

            for (int i = 0; i < n; ++i) {
                vIntBool[i] = vecBool[i];
            }

            cpp11::sexp res = Rf_allocMatrix(LGLSXP, numResults, n);
            int* matBool = LOGICAL(res);

            SerialGlue(matBool, res, vIntBool, nextCmbGrp, nthCmbGrp,
                       nthCmbGrpGmp, FinalTouch, mySample, myBigSamp, startZ,
                       n, numResults, IsArray, IsSample, IsNamed, IsGmp);

            return res;
        } case VecType::Integer : {
            cpp11::sexp res = Rf_allocMatrix(INTSXP, numResults, n);
            int* matInt = INTEGER(res);

            GroupsMain(matInt, res, nextCmbGrp, nthCmbGrp, nthCmbGrpGmp,
                       FinalTouch, vInt, startZ, mySample, myBigSamp, lowerMpz,
                       lower, n, numResults, nThreads, IsArray, IsNamed,
                       Parallel, IsSample, IsGmp);

            if (Rf_isFactor(Rv)) {
                SetFactorClass(res, Rv);
            }

            return res;
        } default : {
            cpp11::sexp res = Rf_allocMatrix(REALSXP, numResults, n);
            double* matNum = REAL(res);

            GroupsMain(matNum, res, nextCmbGrp, nthCmbGrp, nthCmbGrpGmp,
                       FinalTouch, vNum, startZ, mySample, myBigSamp, lowerMpz,
                       lower, n, numResults, nThreads, IsArray, IsNamed,
                       Parallel, IsSample, IsGmp);

            return res;
        }
    }
}
