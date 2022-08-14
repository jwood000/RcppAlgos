#include "Constraints/ConstraintsGeneral.h"
#include "Constraints/ConstraintsSpecial.h"
#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsManager.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/ThreadSafeParts.h"
#include "CombinatoricsResGlue.h"

template <typename T>
void ConstraintsVector(const std::vector<int> &freqs,
                       std::vector<T> &cnstrntVec, std::vector<T> &resVec,
                       std::vector<T> &v, std::vector<T> &tarVals,
                       const std::vector<std::string> &compVec,
                       std::vector<int> &Reps, const std::string &mainFun,
                       const std::string &funTest, std::vector<int> &z,
                       ConstraintType ctype, PartitionType ptype,
                       double lower, mpz_t lowerMpz, int n, int maxRows,
                       int width, int nThreads, bool IsComb, bool IsRep,
                       bool IsMult, bool bUpper, bool xtraCol, bool IsGmp) {

    if (ctype == ConstraintType::General ||
        ctype == ConstraintType::PartitionEsque) {
        ConstraintsGeneral(v, Reps, compVec, cnstrntVec, resVec,
                           tarVals, mainFun, funTest, maxRows, n, width,
                           IsRep, IsComb, IsMult, bUpper, xtraCol, ctype);
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

SEXP ConstraintsReturn(
        const std::vector<int> &freqs, std::vector<double> &vNum,
        std::vector<int> &vInt, std::vector<int> &Reps,
        std::vector<double> &tarVals, std::vector<int> &tarIntVals,
        std::vector<int> &z, const std::vector<std::string> &compVec,
        const std::string &mainFun, const std::string &funTest,
        const PartDesign &part, VecType myType, ConstraintType ctype,
        double userNum, double lower, mpz_t lowerMpz, int n, int m,
        int nRows, int nThreads, double strt, bool IsComb, bool IsRep,
        bool IsMult, bool bUpper, bool xtraCol, bool numUnknown,
        int strtLen, int cap, bool IsGmp
    ) {

    const int width     = (part.isPart) ? part.width : m;
    const int nCols     = (xtraCol) ? width + 1 : width;
    const int lastElem  = n - 1;
    const int lastCol   = width - 1;
    const bool IsVecRet = (numUnknown || part.ptype == PartitionType::Multiset);

    if (part.isPart && !part.solnExist) {
        if (myType == VecType::Integer) {
            cpp11::sexp res = Rf_allocMatrix(INTSXP, 0, width);
            return res;
        } else {
            cpp11::sexp res = Rf_allocMatrix(REALSXP, 0, width);
            return res;
        }
    } else if (myType == VecType::Integer && IsVecRet) {
        std::vector<int> cnstrntVec;
        std::vector<int> resVec;
        const double vecMax = std::floor(cnstrntVec.max_size() / width);
        const double upperBound = std::min(vecMax, dblIntMax);
        const int maxRows = std::min(upperBound, userNum);

        ConstraintsVector(freqs, cnstrntVec, resVec, vInt, tarIntVals,
                          compVec, Reps, mainFun, funTest, z, ctype,
                          part.ptype, lower, lowerMpz, n, maxRows, width,
                          nThreads, IsComb, IsRep, IsMult, bUpper,
                          xtraCol, IsGmp);

        const int vecLen = cnstrntVec.size();
        const int numResult = vecLen / width;

        cpp11::sexp res = Rf_allocMatrix(INTSXP, numResult, nCols);
        int* matInt = INTEGER(res);

        VectorToMatrix(cnstrntVec, resVec, matInt,
                       part.target, numResult, width,
                       upperBound, xtraCol, part.isPart);
        return res;
    } else if (IsVecRet) {
        std::vector<double> cnstrntVec;
        std::vector<double> resVec;
        const double vecMax = std::floor(cnstrntVec.max_size() / width);
        const double upperBound = std::min(vecMax, dblIntMax);
        const int maxRows = std::min(upperBound, userNum);

        ConstraintsVector(freqs, cnstrntVec, resVec, vNum, tarVals,
                          compVec, Reps, mainFun, funTest, z, ctype,
                          part.ptype, lower, lowerMpz, n, maxRows, width,
                          nThreads, IsComb, IsRep, IsMult, bUpper,
                          xtraCol, IsGmp);

        const int vecLen = cnstrntVec.size();
        const int numResult = vecLen / width;

        cpp11::sexp res = Rf_allocMatrix(REALSXP, numResult, nCols);
        double* matDbl = REAL(res);
        VectorToMatrix(cnstrntVec, resVec, matDbl,
                       part.target, numResult, width,
                       upperBound, xtraCol, part.isPart);
        return res;
    } else if (myType == VecType::Integer) {
        cpp11::sexp res = Rf_allocMatrix(INTSXP, nRows, nCols);
        int* matInt = INTEGER(res);

        if (ctype == ConstraintType::PartStandard) {
            StandardPartitions(matInt, z, part.ptype, lower, lowerMpz,
                               nCols, width, nRows, nThreads, lastCol,
                               lastElem, part.mapTar, strtLen, cap, IsRep,
                               IsMult, IsGmp, IsComb, part.includeZero,
                               part.isComp, !part.isWeak);
        } else {
            GeneralPartitions(matInt, vInt, z, part, lower, lowerMpz,
                              nCols, nRows, nThreads, lastCol, lastElem,
                              strtLen, cap, IsComb);
        }

        if (xtraCol) AddResultToParts(matInt, part.target, nRows, width);
        return res;
    } else {
        cpp11::sexp res = Rf_allocMatrix(REALSXP, nRows, nCols);
        double* matDbl = REAL(res);

        GeneralPartitions(matDbl, vNum, z, part, lower, lowerMpz,
                          nCols, nRows, nThreads, lastCol, lastElem,
                          strtLen, cap, IsComb);

        if (xtraCol) AddResultToParts(matDbl, part.target, nRows, width);
        return res;
    }
}

SEXP GetConstraints(
        const PartDesign &part, const std::vector<std::string> &compVec,
        const std::vector<int> &freqs, std::vector<int> &myReps,
        std::vector<double> &vNum, std::vector<int> &vInt,
        std::vector<double> &tarVals, std::vector<int> &tarIntVals,
        std::vector<int> &startZ, const std::string &mainFun,
        const std::string &funTest, funcPtr<double> funDbl, double lower,
        mpz_t lowerMpz, double userNum, ConstraintType ctype, VecType myType,
        int nThreads, int nRows, int n, int strtLen, int cap, int m,
        bool IsComb, bool Parallel, bool IsGmp, bool IsRep, bool IsMult,
        bool bUpper, bool KeepRes, bool numUnknown
    ) {

    if (ctype == ConstraintType::NoConstraint) {
        const int nCols = m + 1;

        if (myType == VecType::Integer) {
            if (!CheckIsInteger(funTest, n, m, vNum, vNum, funDbl,
                                false, IsRep, IsMult, false)) {
                myType = VecType::Numeric;
            }
        }

        if (myType == VecType::Integer) {
            const funcPtr<int> funInt = GetFuncPtr<int>(funTest);
            cpp11::sexp res = Rf_allocMatrix(INTSXP, nRows, nCols);
            int* matInt = INTEGER(res);

            ResultsMain(matInt, vInt, funInt, n, m, IsComb, Parallel,
                        IsRep, IsMult, IsGmp, freqs, startZ, myReps,
                        lower, lowerMpz, nRows, nThreads);

            return res;
        } else {
            funDbl = GetFuncPtr<double>(funTest);
            cpp11::sexp res = Rf_allocMatrix(REALSXP, nRows, nCols);
            double* matNum = REAL(res);

            ResultsMain(matNum, vNum, funDbl, n, m, IsComb, Parallel,
                        IsRep, IsMult, IsGmp, freqs, startZ, myReps,
                        lower, lowerMpz, nRows, nThreads);

            return res;
        }
    } else {
        return ConstraintsReturn(freqs, vNum, vInt, myReps, tarVals,
                                 tarIntVals, startZ, compVec, mainFun,
                                 funTest, part, myType, ctype, userNum,
                                 lower, lowerMpz, n, m, nRows, nThreads,
                                 lower, IsComb, IsRep, IsMult, bUpper,
                                 KeepRes, numUnknown, strtLen, cap, IsGmp);
    }
}
