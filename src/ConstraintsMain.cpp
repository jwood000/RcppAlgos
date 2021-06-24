#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsManager.h"
#include "Partitions/PartitionsUtils.h"
#include "CombinatoricsResGlue.h"
#include "CombinatoricsCnstrt.h"
#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "CheckReturn.h"

template <typename T>
void VectorToMatrix(const std::vector<T> &partsVec, T* mat, int numResult,
                    int width, int partitionLen, int upperBound) {

    for (int count = 0, k = 0; count < numResult; ++count) {
        for (int j = 0; j < width; ++j, ++k) {
            mat[count + numResult * j] = partsVec[k];
        }
    }

    if (partitionLen >= upperBound) {
        Rf_warning("The algorithm terminated early as the number of results "
                       "meeting the criteria exceeds the container's maximum "
                       "capacity or 2^31 - 1");
    }
}

template <typename T>
void ConstraintsVector(const std::vector<T> &v, const std::vector<int> &Reps,
                       std::vector<T> &cnstrntVec, std::vector<int> &z,
                       ConstraintType ctype, PartitionType ptype,
                       int n, int nRows, int width, bool IsComb,
                       bool IsRep, bool IsMult) {

    if (ctype == ConstraintType::General) {

    } else if (ctype == ConstraintType::PartitionEsque) {

    } else if (ctype == ConstraintType::SpecialCnstrnt) {

    } else {
        PartsGenManager(cnstrntVec, v, Reps, z, ptype, width, nRows, IsComb);
    }
}

SEXP ConstraintsReturn(const std::vector<double> &vNum,
                       const std::vector<int> &vInt,
                       const std::vector<int> &Reps, std::vector<int> &z,
                       const PartDesign &part, VecType myType,
                       ConstraintType ctype, double userNumRows, int n,
                       int nRows, int strt, bool IsComb, bool IsRep,
                       bool IsMult) {

    const int width = part.width;

    const auto it = std::find(DistPTypeArr.cbegin(),
                              DistPTypeArr.cend(), part.ptype);

    const bool partsUnknown = (IsMult && it == DistPTypeArr.cend()) ||
                (part.ptype == PartitionType::RepCapped && !IsComb);

    const bool numUnknown = ctype == ConstraintType::PartitionEsque ||
                            ctype == ConstraintType::SpecialCnstrnt ||
                            ctype == ConstraintType::General        ||
                            partsUnknown;

    if (numUnknown) {
        constexpr double int_max = std::numeric_limits<int>::max();
        if (myType == VecType::Integer) {
            std::vector<int> cnstrntVec;
            const double vecMax = cnstrntVec.max_size() / width;
            const double upperBound = std::min(vecMax, int_max);
            const int maxRows = std::min(upperBound, userNumRows);

            ConstraintsVector(vInt, Reps, cnstrntVec, z, ctype,
                              part.ptype, n, maxRows,
                              width, IsComb, IsRep, IsMult);

            const int vecLen = cnstrntVec.size();
            const int numResult = vecLen / width;

            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, numResult, width));
            int* matInt = INTEGER(res);
            VectorToMatrix(cnstrntVec, matInt, numResult,
                           width, vecLen, upperBound);
            UNPROTECT(1);
            return res;
        } else {
            std::vector<double> cnstrntVec;
            const double vecMax = cnstrntVec.max_size() / width;
            const double upperBound = std::min(vecMax, int_max);
            const int maxRows = std::min(upperBound, userNumRows);

            ConstraintsVector(vNum, Reps, cnstrntVec, z, ctype,
                              part.ptype, n, maxRows,
                              width, IsComb, IsRep, IsMult);

            const int vecLen = cnstrntVec.size();
            const int numResult = vecLen / width;

            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, numResult, width));
            double* matDbl = REAL(res);
            VectorToMatrix(cnstrntVec, matDbl, numResult,
                           width, vecLen, upperBound);
            UNPROTECT(1);
            return res;
        }
    } else {
        const int lastCol = width - 1;
        const int lastElem = n - 1;
        int boundary = lastCol;
        int edge = boundary - 1;

        if (myType == VecType::Integer) {
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, nRows, width));
            int* matInt = INTEGER(res);

            if (ctype == ConstraintType::PartStandard) {
                PartsStdManager(matInt, z, width, lastCol,
                                boundary, edge, strt, nRows, IsComb, IsRep);
            } else {
                PartsGenManager(matInt, part, vInt, z, width,
                                lastCol, lastElem, boundary, edge,
                                strt, nRows, IsRep, IsComb);
            }

            UNPROTECT(1);
            return res;
        } else {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, nRows, width));
            double* matDbl = REAL(res);

            PartsGenManager(matDbl, part, vNum, z, width,
                            lastCol, lastElem, boundary, edge,
                            strt, nRows, IsRep, IsComb);

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

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> targetIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compFunVec;
    std::vector<double> targetVals;

    ConstraintType ctype = ConstraintType::NoConstraint;
    PartDesign part;

    part.isRep = IsRep;
    part.isMult = IsMult;
    part.mIsNull = Rf_isNull(Rm);

    if (IsConstrained) {
        ConstraintSetup(vNum, myReps, targetVals, targetIntVals, funDbl,
                        part, ctype, n, m, compFunVec, mainFun, myType,
                        Rtarget, RcompFun, Rtolerance, Rlow, IsComb, false);
    }

    const double computedRows = (part.count > 0) ? part.count :
        GetComputedRows(IsMult, IsComb, IsRep, n, m, Rm, freqs, myReps);

    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp && part.count == 0) {
        GetComputedRowMpz(computedRowsMpz, IsMult,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }

    double lower = 0;
    double upper = 0;

    bool bLower = false;
    bool bUpper = false;

    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);

    mpz_init(lowerMpz[0]);
    mpz_init(upperMpz[0]);

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz.get(), upperMpz.get(), computedRowsMpz, computedRows);

    std::vector<int> startZ(m);

    if (ctype < ConstraintType::PartMapping) {
        SetStartZ(myReps, freqs, startZ, IsComb, n, m,
                  lower, lowerMpz[0], IsRep, IsMult, IsGmp);
    } else {
        startZ = part.startZ;
    }

    // This is used when we are unable to calculate the number of results
    // upfront (E.g. comboGeneral(rnorm(10), 5, constraintFun = "sum,
    //                            comparisonFun = "<=", limitConstraints = 1))
    double userNumRows = 0;

    // This variable is used in determining the number of results. If the
    // output is constrained and the ConstraintType is "General" or
    // "PartitionEsque", it means we really don't know how many results
    // we have. The computedRows above is a strict upper bound but not
    // necessarily the least upper bound. In these cases, we don't want
    // to unnecessarily throw an error when computedRows exceeds 2^31 - 1.
    const bool IsGenCnstrd = IsConstrained &&
                             (ctype == ConstraintType::General ||
                              ctype == ConstraintType::PartitionEsque ||
                              part.ptype == PartitionType::Multiset);

    SetNumResults(IsGmp, bLower, bUpper, IsGenCnstrd, upperMpz.get(),
                  lowerMpz.get(), lower, upper, computedRows,
                  computedRowsMpz, nRows, userNumRows);

    if (ctype == ConstraintType::NoConstraint) {
        int nThreads = 1;
        int maxThreads = 1;
        CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                       VecType::Integer, "maxThreads");

        int nCols = m + 1;
        const int limit = 20000;
        SetThreads(Parallel, maxThreads, nRows, myType, nThreads, RnThreads, limit);

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
        return ConstraintsReturn(vNum, vInt, myReps, startZ, part,
                                 myType, ctype, userNumRows, n, nRows,
                                 lower, IsComb, IsRep, IsMult);
    }
}
