#include "Constraints/PartitionsEsqueAlgo.h"
#include "Constraints/ConstraintsGeneral.h"
#include "Constraints/ConstraintsSpecial.h"
#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsManager.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/NthPartition.h"
#include "CombinatoricsResGlue.h"
#include "CombinatoricsCnstrt.h"
#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "CheckReturn.h"

template <typename T>
void AddResultToParts(T* mat, std::int64_t result,
                      std::size_t numResult,
                      std::size_t width) {

    const T t_result = result;
    const std::size_t limit = static_cast<std::size_t>(numResult) *
                              static_cast<std::size_t>(width + 1);

    for (std::size_t i = numResult * width; i < limit; ++i) {
        mat[i] = t_result;
    }
}

template <typename T>
void VectorToMatrix(const std::vector<T> &cnstrntVec,
                    const std::vector<T> &resVec, T* mat,
                    std::int64_t result, std::size_t numResult,
                    std::size_t width, int upperBound,
                    bool xtraCol, bool IsPart) {

    if (numResult >= (upperBound - 1)) {
        Rf_warning("The algorithm terminated early as the number of results"
                   " meeting the criteria exceeds the container's maximum"
                   " capacity or 2^31 - 1");
    }

    for (std::size_t count = 0, k = 0; count < numResult; ++count) {
        for (std::size_t j = 0; j < width; ++j, ++k) {
            mat[count + numResult * j] = cnstrntVec[k];
        }
    }

    if (xtraCol) {
        const std::size_t limit = static_cast<std::size_t>(numResult) *
                                  static_cast<std::size_t>(width + 1);

        if (IsPart) {
            AddResultToParts(mat, result, numResult, width);
        } else {
            for (std::size_t i = numResult * width, k = 0;
                 i < limit; ++i, ++k) {

                mat[i] = resVec[k];
            }
        }
    }
}

template <typename T>
void ConstraintsVector(const std::vector<int> &freqs,
                       std::vector<T> &cnstrntVec, std::vector<T> &resVec,
                       std::vector<T> &v, std::vector<T> &tarVals,
                       const std::vector<std::string> &compVec,
                       std::vector<int> &Reps, const std::string &mainFun,
                       std::vector<int> &z, ConstraintType ctype,
                       PartitionType ptype, double lower, mpz_t lowerMpz,
                       int n, int maxRows, int width, int nThreads,
                       bool IsComb, bool IsRep, bool IsMult,
                       bool bUpper, bool xtraCol, bool IsGmp) {

    if (ctype == ConstraintType::General) {
        ConstraintsGeneral(v, Reps, compVec, cnstrntVec, resVec,
                           tarVals, mainFun, maxRows, n, width, IsRep,
                           IsComb, IsMult, bUpper, xtraCol);
    } else if (ctype == ConstraintType::PartitionEsque) {
        PartitionsEsqueAlgo(v, tarVals, Reps, mainFun, compVec.front(),
                            cnstrntVec, resVec, maxRows, n, width, IsRep,
                            IsComb, xtraCol, IsMult, bUpper);
    } else if (ctype == ConstraintType::SpecialCnstrnt) {
        ConstraintsSpecial(v, tarVals, compVec, Reps, freqs, cnstrntVec,
                           resVec, mainFun, z, lower, lowerMpz, n, width,
                           maxRows, nThreads, IsRep, xtraCol, IsComb,
                           IsMult, IsGmp);
    } else {
        PartsGenManager(cnstrntVec, v, Reps, z,
                        ptype, width, maxRows, IsComb);
    }
}

SEXP ConstraintsReturn(const std::vector<int> &freqs,
                       std::vector<double> &vNum, std::vector<int> &vInt,
                       std::vector<int> &Reps, std::vector<double> &tarVals,
                       std::vector<int> &tarIntVals, std::vector<int> &z,
                       const std::vector<std::string> &compVec,
                       const std::string &mainFun, const PartDesign &part,
                       VecType myType, ConstraintType ctype, double userNum,
                       double lower, mpz_t lowerMpz, int n, int m, int nRows,
                       int nThreads, double strt, bool IsComb, bool IsRep,
                       bool IsMult, bool bUpper, bool xtraCol,
                       bool numUnknown, bool IsGmp) {

    const int width = (part.isPart) ? part.width : m;
    const int nCols = (xtraCol) ? width + 1 : width;

    if (part.isPart && !part.solnExist) {
        if (myType == VecType::Integer) {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, 0, width));
            UNPROTECT(1);
            return res;
        } else {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, 0, width));
            UNPROTECT(1);
            return res;
        }
    } else if (numUnknown) {
        if (myType == VecType::Integer) {
            std::vector<int> cnstrntVec;
            std::vector<int> resVec;
            const double vecMax = std::floor(cnstrntVec.max_size() / width);
            const double upperBound = std::min(vecMax, dblIntMax);
            const int maxRows = std::min(upperBound, userNum);

            ConstraintsVector(freqs, cnstrntVec, resVec, vInt, tarIntVals,
                              compVec, Reps, mainFun, z, ctype, part.ptype,
                              lower, lowerMpz, n, maxRows, width, nThreads,
                              IsComb, IsRep, IsMult, bUpper, xtraCol, IsGmp);

            const int vecLen = cnstrntVec.size();
            const int numResult = vecLen / width;

            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, numResult, nCols));
            int* matInt = INTEGER(res);
            VectorToMatrix(cnstrntVec, resVec, matInt,
                           part.target, numResult, width,
                           upperBound, xtraCol, part.isPart);
            UNPROTECT(1);
            return res;
        } else {
            std::vector<double> cnstrntVec;
            std::vector<double> resVec;
            const double vecMax = std::floor(cnstrntVec.max_size() / width);
            const double upperBound = std::min(vecMax, dblIntMax);
            const int maxRows = std::min(upperBound, userNum);

            ConstraintsVector(freqs, cnstrntVec, resVec, vNum, tarVals,
                              compVec, Reps, mainFun, z, ctype, part.ptype,
                              lower, lowerMpz, n, maxRows, width, nThreads,
                              IsComb, IsRep, IsMult, bUpper, xtraCol, IsGmp);

            const int vecLen = cnstrntVec.size();
            const int numResult = vecLen / width;

            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, numResult, nCols));
            double* matDbl = REAL(res);
            VectorToMatrix(cnstrntVec, resVec, matDbl,
                           part.target, numResult, width,
                           upperBound, xtraCol, part.isPart);
            UNPROTECT(1);
            return res;
        }
    } else {
        const int lastElem = n - 1;
        const int lastCol = width - 1;

        if (myType == VecType::Integer) {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, nRows, nCols));
            int* matInt = INTEGER(res);

            if (ctype == ConstraintType::PartStandard) {
                PartsStdManager(matInt, z, width, lastElem,
                                lastCol, nRows, IsComb, IsRep);
            } else {
                PartsGenManager(matInt, part, vInt, z, width, lastElem,
                                lastCol, nRows, IsComb, IsRep);
            }

            if (xtraCol) AddResultToParts(matInt, part.target, nRows, width);
            UNPROTECT(1);
            return res;
        } else {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, nRows, nCols));
            double* matDbl = REAL(res);

            PartsGenManager(matDbl, part, vNum, z, width, lastElem,
                            lastCol, nRows, IsComb, IsRep);

            if (xtraCol) AddResultToParts(matDbl, part.target, nRows, width);
            UNPROTECT(1);
            return res;
        }
    }
}

SEXP CombinatoricsCnstrt(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs,
                         SEXP Rlow, SEXP Rhigh, SEXP RmainFun,
                         SEXP RcompFun, SEXP Rtarget, SEXP RIsComb,
                         SEXP RKeepRes, SEXP Rparallel, SEXP RnThreads,
                         SEXP RmaxThreads, SEXP Rtolerance) {

    int n = 0;
    int m = 0;
    int nRows = 0;

    bool IsMult = false;
    VecType myType = VecType::Integer;

    std::vector<double> vNum;
    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;

    bool KeepRes  = CleanConvert::convertLogical(RKeepRes, "keepResults");
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRep    = CleanConvert::convertLogical(RisRep, "repetition");

    const bool IsComb = CleanConvert::convertLogical(RIsComb, "IsComb");
    const bool IsConstrained = CheckConstrnd(RmainFun, RcompFun, Rtarget);

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum, Rv,
              RFreqs, Rm, n, m, IsMult, IsRep, IsConstrained);

    if (!Rf_isString(RmainFun) || Rf_length(RmainFun) != 1) {
        Rf_error("contraintFun must be one of the following:"
                 " 'prod', 'sum', 'mean', 'max', or 'min'");
    }

    const std::string mainFun(CHAR(STRING_ELT(RmainFun, 0)));
    const auto funIt = std::find(mainFunSet.begin(), mainFunSet.end(), mainFun);

    if (funIt == mainFunSet.end()) {
        Rf_error("contraintFun must be one of the following:"
                 " 'prod', 'sum', 'mean', 'max', or 'min'");
    }

    // Must be defined inside IsInteger check as tarVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> tarIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compVec;
    std::vector<double> tarVals;

    ConstraintType ctype = ConstraintType::NoConstraint;
    PartDesign part;

    part.isRep = IsRep;
    part.isMult = IsMult;
    part.mIsNull = Rf_isNull(Rm);

    if (IsConstrained) {
        ConstraintSetup(vNum, myReps, tarVals, tarIntVals, funDbl,
                        part, ctype, n, m, compVec, mainFun, myType,
                        Rtarget, RcompFun, Rtolerance, Rlow, IsComb, false);
    }

    const double computedRows = (part.count > 0) ? part.count :
        GetComputedRows(IsMult, IsComb, IsRep, n, m, Rm, freqs, myReps);

    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp && part.isPart) {
        mpz_set(computedRowsMpz, part.bigCount);
    } else if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }

    // This variable is used in determining the number of results. If the
    // output is constrained and the ConstraintType is "General" or
    // "PartitionEsque", it means we really don't know how many results
    // we have. The computedRows above is a strict upper bound but not
    // necessarily the least upper bound. In these cases, we don't want
    // to unnecessarily throw an error when computedRows exceeds 2^31 - 1.
    const bool numUnknown = ctype == ConstraintType::PartitionEsque ||
                            ctype == ConstraintType::SpecialCnstrnt ||
                            ctype == ConstraintType::General        ||
                            part.numUnknown;

    double lower = 0;
    double upper = 0;

    bool bLower = false;
    bool bUpper = false;

    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);

    mpz_init(lowerMpz[0]);
    mpz_init(upperMpz[0]);

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz.get(), upperMpz.get(), computedRowsMpz,
              computedRows);

    std::vector<int> startZ(m);

    if (ctype < ConstraintType::PartMapping) {
        SetStartZ(myReps, freqs, startZ, IsComb, n, m,
                  lower, lowerMpz[0], IsRep, IsMult, IsGmp);
    } else {
        if (bLower) {
            const int strtLen = std::count_if(part.startZ.cbegin(),
                                              part.startZ.cend(),
                                              [](int i){return i > 0;});

            const int k = (part.ptype == PartitionType::DstctSpecial ||
                           part.ptype == PartitionType::DstctStdAll) ?
                           strtLen : n;

            const nthPartsPtr nthPartFun = GetNthPartsFunc(part.ptype,
                                                           IsGmp);
            startZ = nthPartFun(part.mapTar, part.width,
                                k, lower, lowerMpz[0]);

            if (ctype == ConstraintType::PartStandard && !part.includeZero) {
                for (auto &z_i: startZ) {
                    ++z_i;
                }
            }
        } else {
            startZ = part.startZ;
        }
    }

    // This is used when we are unable to calculate the number of results
    // upfront (E.g. comboGeneral(rnorm(10), 5, constraintFun = "sum,
    //                            comparisonFun = "<=", limitConstraints = 1))
    double userNum = 0;
    const bool bSetNum = !numUnknown ||
        ctype == ConstraintType::SpecialCnstrnt;

    SetNumResults(IsGmp, bLower, bUpper, bSetNum, upperMpz.get(),
                  lowerMpz.get(), lower, upper, computedRows,
                  computedRowsMpz, nRows, userNum);

    int nThreads = 1;
    int maxThreads = 1;
    CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                   VecType::Integer, "maxThreads");

    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows,
               myType, nThreads, RnThreads, limit);

    if (ctype == ConstraintType::NoConstraint) {
        const int nCols = m + 1;

        if (myType == VecType::Integer) {
            if (!CheckIsInteger(mainFun, n, m, vNum, vNum, funDbl)) {
                myType = VecType::Numeric;
            }
        }

        if (myType == VecType::Integer) {
            const funcPtr<int> funInt = GetFuncPtr<int>(mainFun);
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, nRows, nCols));
            int* matInt = INTEGER(res);

            ResultsMain(matInt, vInt, funInt, n, m, IsComb, Parallel,
                        IsRep, IsMult, IsGmp, freqs, startZ, myReps,
                        lower, lowerMpz[0], nRows, nThreads);

            UNPROTECT(1);
            return res;
        } else {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, nRows, nCols));
            double* matNum = REAL(res);

            ResultsMain(matNum, vNum, funDbl, n, m, IsComb, Parallel,
                        IsRep, IsMult, IsGmp, freqs, startZ, myReps,
                        lower, lowerMpz[0], nRows, nThreads);

            UNPROTECT(1);
            return res;
        }
    } else {
        return ConstraintsReturn(freqs, vNum, vInt, myReps, tarVals,
                                 tarIntVals, startZ, compVec, mainFun,
                                 part, myType, ctype, userNum, lower,
                                 lowerMpz[0], n, m, nRows, nThreads, lower,
                                 IsComb, IsRep, IsMult, bUpper, KeepRes,
                                 numUnknown, IsGmp);
    }
}
