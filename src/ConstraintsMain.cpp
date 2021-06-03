#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsUtils.h"
#include "CombinatoricsResGlue.h"
#include "CombinatoricsCnstrt.h"
#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "CheckReturn.h"

// template <typename typeRcpp, typename typeVector>
// typeRcpp ConstraintReturn(int n, int m, const std::string &mainFun, const std::vector<std::string> &compFunVec,
//                           const std::vector<typeVector> &targetVals, std::vector<typeVector> &v, bool bLower,
//                           double lower, bool bUserRows, double computedRows, bool IsRep, int nRows, bool KeepRes,
//                           std::vector<int> &z, bool IsMult, bool IsComb, bool mIsNull, double userNumRows,
//                           std::vector<int> &myReps, const std::vector<int> &freqs, PartitionType PartType,
//                           ConstraintType ConstType, bool distGetAll) {
//     
//     if (SpecialCase) {
//         return ConstraintsSpecial<typeRcpp>(n, m, v, IsRep, nRows, KeepRes, z, lower, mainFun, 
//                                             IsMult, computedRows, compFunVec, targetVals, IsComb,
//                                             freqs, bLower, userNumRows);
//     }
//     
//     // For bool bUserRows, we pass bUpper as we know bLower must be false (See CheckSpecialCase)
//     if (ConstType > ConstraintType::PartitionEsque) {
//         return Partitions::PartitionsMain<typeRcpp>(v, z, myReps, PartType, targetVals[0], n, m, IsRep,
//                                                     IsMult, userNumRows, IsComb, KeepRes, bUserRows,
//                                                     mIsNull, distGetAll);
//     } else if (ConstType == ConstraintType::PartitionEsque) {
//         return PartitionEsqueAlgo<typeRcpp>(n, m, v, IsRep, mainFun, compFunVec.front(), targetVals,
//                                             userNumRows, IsComb, KeepRes, myReps, IsMult, bUserRows);
//     } else {
//         return CombinatoricsConstraints<typeRcpp>(n, m, v, IsRep, mainFun, compFunVec, targetVals,
//                                                   userNumRows, IsComb, KeepRes, myReps, IsMult, bUserRows);
//     }
// }

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
    
    ConstraintType ctype;
    PartDesign part;

    part.isRep = IsRep;
    part.isMult = IsMult;
    part.mIsNull = Rf_isNull(Rm);
    
    ConstraintSetup(vNum, myReps, targetVals, targetIntVals,
                    funDbl, part, ctype, n, m, compFunVec,
                    mainFun, myType, Rtarget, RcompFun, Rtolerance,
                    Rlow, IsConstrained, false);
    
    
    // const bool IsPartition = CheckPartition(compFunVec, vNum, mainFun,
    //                                         targetVals, part, Rlow, lenV,
    //                                         m, tolerance, IsBetweenComp);
    // 
    // if (IsPartition) {
    //     GetPartitionDesign(Reps, vNum, part, ctype, lenV, m, bCalcMulti);
    // }
    
    if (ctype >= ConstraintType::PartMapping) {
        
    } else {
        
    }
    
    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);

    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp) {
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
    
    if (ctype == ConstraintType::General) {
        SetStartZ(myReps, freqs, startZ, IsComb, n, m,
                  lower, lowerMpz[0], IsRep, IsMult, IsGmp);
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
                              ctype == ConstraintType::PartitionEsque);
    
    SetNumResults(IsGmp, bLower, bUpper, IsGenCnstrd, upperMpz.get(),
                  lowerMpz.get(), lower, upper, computedRows,
                  computedRowsMpz, nRows, userNumRows);

    // if (IsConstrained) {
    //     if (myType == VecType::Integer) {
    //         return ConstraintReturn<Rcpp::IntegerMatrix>(n, m, mainFun, compFunVec, targetIntVals, vInt, bLower,
    //                                                      lower, bUpper, computedRows, IsRep, nRows, KeepRes, startZ,
    //                                                      IsMult, IsComb, Rf_isNull(Rm), userNumRows, myReps, freqs,
    //                                                      PartType, ConstType, distinctTest.getAll);
    //     } else {
    //         return ConstraintReturn<Rcpp::NumericMatrix>(n, m, mainFun, compFunVec, targetVals, vNum, bLower,
    //                                                      lower, bUpper, computedRows, IsRep, nRows, KeepRes, startZ,
    //                                                      IsMult, IsComb, Rf_isNull(Rm), userNumRows, myReps, freqs,
    //                                                      PartType, ConstType, distinctTest.getAll);
    //     }
    // } else {
        int nThreads = 1;
        int maxThreads = 1;
        CleanConvert::convertPrimitive(RmaxThreads, maxThreads,
                                       VecType::Integer, "maxThreads");

        int nCols = m + 1;
        const int limit = 20000;
        SetThreads(Parallel, maxThreads, nRows, myType, nThreads, RnThreads, limit);

        if (myType == VecType::Integer)
            if (!CheckIsInteger(mainFun, n, m, vNum, vNum, funDbl))
                myType = VecType::Numeric;

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
    // }
}
