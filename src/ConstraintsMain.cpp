#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionEnums.h"
#include "CombinatoricsResGlue.h"
#include "CombinatoricsCnstrt.h"
#include "Cpp14MakeUnique.h"
#include "ComputedCount.h"
#include "CheckReturn.h"
#include "SetUpUtils.h"

// #include "CombPermResultPtr.h"
// #include "ConstraintsGeneral.h"
// #include "ConstraintsUtils.h"
// #include "PartitionsMain.h"
// #include "PartitionEsqueAlgo.h"
// #include "ConstraintsSpecial.h"
// #include "PartitionsCounts.h"
// #include "CheckStdRet.h"
// #include "RMatrix.h"
// #include <RcppThread/ThreadPool.hpp>
// 
// template <typename typeRcpp, typename typeVector>
// typeRcpp ConstraintReturn(int n, int m, const std::string &mainFun, const std::vector<std::string> &compFunVec,
//                           const std::vector<typeVector> &targetVals, std::vector<typeVector> &v, bool bLower,
//                           double lower, bool bUserRows, double computedRows, bool IsRep, int nRows, bool KeepRes,
//                           std::vector<int> &z, bool IsMult, bool IsComb, bool mIsNull, double userNumRows,
//                           std::vector<int> &myReps, const std::vector<int> &freqs, PartitionType PartType,
//                           ConstraintType ConstType, bool distGetAll) {
//     
//     const bool SpecialCase = CheckSpecialCase(n, bLower, mainFun, v);
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

SEXP CombinatoricsCnstrt(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow,
                         SEXP Rhigh, SEXP f1, SEXP f2, SEXP Rtarget, SEXP RIsComb,
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

    bool KeepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    const bool IsComb = CleanConvert::convertLogical(RIsComb, "IsComb");

    SetType(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    const bool IsConstrained = CheckConstrnd(f1, f2, Rtarget);

    if (IsConstrained) {
        for (int i = (vNum.size() - 1); i >= 0; --i)
            if (ISNAN(vNum[i]))
                vNum.erase(vNum.begin() + i);

        n = vNum.size();
    }

    SetFreqsAndM(RFreqs, IsMult, myReps, IsRep, freqs, Rm, n, m);
    const char* f1_chars;
    f1_chars = CHAR(STRING_ELT(f1, 0));
    
    const std::string mainFun(f1_chars);
    const auto funIt = std::find(mainFunSet.begin(), mainFunSet.end(), mainFun);

    if (funIt == mainFunSet.end()) {
        Rf_error("contraintFun must be one of the following: prod, sum, mean, max, or min");
    }

    std::vector<std::string> compFunVec;
    std::vector<double> targetVals;

    ConstraintType ConstType = ConstraintType::General;
    PartitionType PartType = PartitionType::Traditional;
    // distinctType distinctTest;

    const bool allNeg = std::all_of(vNum.cbegin(), vNum.cend(),
                                    [](double v_i) {return v_i <= 0;});
    const bool allPos = std::all_of(vNum.cbegin(), vNum.cend(),
                                    [](double v_i) {return v_i >= 0;});
    const Sign mySign = allPos ? Sign::Positive :
                       (allNeg ? Sign::Negitive : Sign::MixedBag);

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> targetIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    // if (IsConstrained) {                // numOnly = true, checkWhole = false, negPoss = true
    //     CleanConvert::convertVector(Rtarget, targetVals, "limitConstraints", true, false, true);
    //     compFunVec = Rcpp::as<std::vector<std::string>>(f2);
    // 
    //     bool IsBetweenComp = false;
    //     ConstraintSetup(compFunVec, targetVals, IsBetweenComp);
    // 
    //     if (myType == VecType::Integer && !CheckIsInteger(mainFun, n, m, vNum, targetVals, funDbl, true)) {
    //         myType = VecType::Numeric;
    //     }
    // 
    //     double tolerance = 0;
    //     AdjustTargetVals(n, myType, targetVals, targetIntVals,
    //                      Rtolerance, compFunVec, tolerance, mainFun, vNum);
    // 
    //     if (myType == VecType::Integer) {
    //         GetPartitionCase(compFunVec, vInt, mainFun, targetIntVals, mySign,
    //                          PartType, ConstType, distinctTest, Rlow, myReps, n, m,
    //                          tolerance, IsMult, IsRep, IsBetweenComp, Rf_isNull(Rm));
    //     } else {
    //         GetPartitionCase(compFunVec, vNum, mainFun, targetVals, mySign,
    //                          PartType, ConstType, distinctTest, Rlow, myReps, n, m,
    //                          tolerance, IsMult, IsRep, IsBetweenComp, Rf_isNull(Rm));
    //     }
    // }

    std::vector<int> startZ(m);
    double computedRows = 0;

    // Rcpp::Rcout << static_cast<std::underlying_type<PartitionType>::type>(PartType) << std::endl;
    // Rcpp::Rcout << static_cast<std::underlying_type<ConstraintType>::type>(ConstType) << std::endl;
    // 
    // switch (PartType) {
    //     case PartitionType::Traditional : {
    //         Rcpp::print(Rcpp::wrap("Traditional"));
    //         break;
    //     } case PartitionType::TradNoZero : {
    //         Rcpp::print(Rcpp::wrap("TradNoZero"));
    //         break;
    //     } case PartitionType::TradCapped : {
    //         Rcpp::print(Rcpp::wrap("TradCapped"));
    //         break;
    //     } case PartitionType::DstctStdAll : {
    //         Rcpp::print(Rcpp::wrap("DstctStdAll"));
    //         break;
    //     } case PartitionType::DstctShort : {
    //         Rcpp::print(Rcpp::wrap("DstctShort"));
    //         break;
    //     } case PartitionType::DstctSpecial : {
    //         Rcpp::print(Rcpp::wrap("DstctSpecial"));
    //         break;
    //     } case PartitionType::DstctOneZero : {
    //         Rcpp::print(Rcpp::wrap("DstctOneZero"));
    //         break;
    //     } case PartitionType::DstctNoZero : {
    //         Rcpp::print(Rcpp::wrap("DstctNoZero"));
    //         break;
    //     } case PartitionType::DistCapped : {
    //         Rcpp::print(Rcpp::wrap("DistCapped"));
    //         break;
    //     }
    // }
    // 
    // if (ConstType > ConstraintType::General) {
    //     // vNum and myReps were sorted in GetPartitionCase
    //     bool IncludeZero = (vNum.front() == 0);
    //     int targetInt = static_cast<int>(targetVals[0]);
    //     const std::vector<std::int64_t> v64(vNum.cbegin(), vNum.cend());
    // 
    //     SetStartPartitionZ(PartType, distinctTest, v64, startZ, myReps,
    //                        targetInt, n, m, IncludeZero, IsRep, IsMult);
    //     return Rcpp::wrap(startZ);
    //     if (startZ.empty()) {
    //         if (myType == VecType::Integer) {
    //             Rcpp::IntegerMatrix trivialIntRet(0, m);
    //             return trivialIntRet;
    //         } else {
    //             Rcpp::NumericMatrix trivialNumRet(0, m);
    //             return trivialNumRet;
    //         }
    //     }
    // 
    //     computedRows = GetComputedPartsComps(startZ, PartType, targetInt, m,
    //                                          IsComb, IncludeZero, Rf_isNull(Rm));
    // } else {
        computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                       n, m, Rm, freqs, myReps);
    // }
    // 
    // return Rcpp::wrap(ConstType > ConstraintType::General);

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

    // if (ConstType == ConstraintType::General) {
        SetStartZ(myReps, freqs, startZ, IsComb, n, m,
                  lower, lowerMpz[0], IsRep, IsMult, IsGmp);
    // }

    double userNumRows = 0;
    const bool IsGenCnstrd = false; //(IsConstrained && ConstType == ConstraintType::General);
    
    SetNumResults(IsGmp, bLower, bUpper, IsGenCnstrd, upperMpz.get(), lowerMpz.get(),
                  lower, upper, computedRows, computedRowsMpz, nRows, userNumRows);

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
        CleanConvert::convertPrimitive(RmaxThreads, maxThreads, "maxThreads");
        
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
