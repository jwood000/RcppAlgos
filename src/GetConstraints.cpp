#include "Constraints/PartitionsEsqueAlgo.h"
#include "Constraints/ConstraintsGeneral.h"
#include "Constraints/ConstraintsSpecial.h"
#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsManager.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/ThreadSafeParts.h"
#include "CombinatoricsResGlue.h"

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

SEXP ConstraintsReturn(
        const std::vector<int> &freqs, std::vector<double> &vNum,
        std::vector<int> &vInt, std::vector<int> &Reps,
        std::vector<double> &tarVals, std::vector<int> &tarIntVals,
        std::vector<int> &z, const std::vector<std::string> &compVec,
        const std::string &mainFun, const PartDesign &part, VecType myType,
        ConstraintType ctype, double userNum, double lower, mpz_t lowerMpz,
        int n, int m, int nRows, int nThreads, double strt, bool IsComb,
        bool IsRep, bool IsMult, bool bUpper, bool xtraCol, bool numUnknown,
        int strtLen, int cap, bool IsGmp
) {

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
    } else if (numUnknown || part.ptype == PartitionType::Multiset) {
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
                StandardPartitions(matInt, z, part.ptype, lower, lowerMpz,
                                   nCols, width, nRows, nThreads, lastCol,
                                   lastElem, part.mapTar, strtLen, cap, IsRep,
                                   IsMult, IsGmp, IsComb, part.includeZero);
            } else {
                GeneralPartitions(matInt, vInt, z, part, lower, lowerMpz,
                                  nCols, nRows, nThreads, lastCol, lastElem,
                                  strtLen, cap, IsComb);
            }

            if (xtraCol) AddResultToParts(matInt, part.target, nRows, width);
            UNPROTECT(1);
            return res;
        } else {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, nRows, nCols));
            double* matDbl = REAL(res);

            GeneralPartitions(matDbl, vNum, z, part, lower, lowerMpz,
                              nCols, nRows, nThreads, lastCol, lastElem,
                              strtLen, cap, IsComb);

            if (xtraCol) AddResultToParts(matDbl, part.target, nRows, width);
            UNPROTECT(1);
            return res;
        }
    }
}

SEXP GetConstraints(
        const PartDesign &part, const std::vector<std::string> &compVec,
        const std::vector<int> &freqs, std::vector<int> &myReps,
        std::vector<double> &vNum, std::vector<int> &vInt,
        std::vector<double> &tarVals, std::vector<int> &tarIntVals,
        std::vector<int> &startZ, const std::string &mainFun,
        funcPtr<double> funDbl, double lower, mpz_t lowerMpz, double userNum,
        ConstraintType ctype, VecType myType, int nThreads, int nRows, int n,
        int strtLen, int cap, int m, bool IsComb, bool Parallel, bool IsGmp,
        bool IsRep, bool IsMult, bool bUpper, bool KeepRes, bool numUnknown
    ) {

    if (ctype == ConstraintType::NoConstraint) {
        const int nCols = m + 1;

        if (myType == VecType::Integer) {
            if (!CheckIsInteger(mainFun, n, m, vNum, vNum, funDbl,
                                false, IsRep, IsMult, false)) {
                myType = VecType::Numeric;
            }
        }

        if (myType == VecType::Integer) {
            const funcPtr<int> funInt = GetFuncPtr<int>(mainFun);
            SEXP res = PROTECT(Rf_allocMatrix(INTSXP, nRows, nCols));
            int* matInt = INTEGER(res);

            ResultsMain(matInt, vInt, funInt, n, m, IsComb, Parallel,
                        IsRep, IsMult, IsGmp, freqs, startZ, myReps,
                        lower, lowerMpz, nRows, nThreads);

            UNPROTECT(1);
            return res;
        } else {
            SEXP res = PROTECT(Rf_allocMatrix(REALSXP, nRows, nCols));
            double* matNum = REAL(res);

            ResultsMain(matNum, vNum, funDbl, n, m, IsComb, Parallel,
                        IsRep, IsMult, IsGmp, freqs, startZ, myReps,
                        lower, lowerMpz, nRows, nThreads);

            UNPROTECT(1);
            return res;
        }
    } else {
        return ConstraintsReturn(freqs, vNum, vInt, myReps, tarVals,
                                 tarIntVals, startZ, compVec, mainFun,
                                 part, myType, ctype, userNum, lower,
                                 lowerMpz, n, m, nRows, nThreads, lower,
                                 IsComb, IsRep, IsMult, bUpper, KeepRes,
                                 numUnknown, strtLen, cap, IsGmp);
    }
}
