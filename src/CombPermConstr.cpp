#include "ConstraintsMaster.h"
#include "GmpCombPermUtils.h"
#include "PartitionsCounts.h"
#include "RMatrix.h"
#include <RcppThread.h>

template <typename typeRcpp, typename T>
void MasterResRet(typeRcpp &matRcpp, const std::vector<T> &v, funcPtr<T> myFun, int n, int m,
                  bool IsRep, bool IsComb, bool IsMult, bool IsGmp, const std::vector<int> &freqs,
                  std::vector<int> z, const std::vector<int> &myReps, double lower, mpz_t &lowerMpz,
                  int nRows, int nThreads, bool Parallel) {

    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(matRcpp);
        RcppThread::ThreadPool pool(nThreads);
        const int stepSize = nRows / nThreads;
        int nextStep = stepSize;
        int step = 0;
        
        Rcpp::XPtr<combPermResPtr<RcppParallel::RMatrix<T>, T>>
            xpFunCoPeResPtr = putCombResPtrInXPtr<RcppParallel::RMatrix<T>, T>(IsComb, IsMult, IsRep);
        const combPermResPtr<RcppParallel::RMatrix<T>, T> myFunParCombPerm = *xpFunCoPeResPtr;
        
        Rcpp::XPtr<nthResultPtr> xpNthComb = putNthResPtrInXPtr(IsComb, IsMult, IsRep, IsGmp);
        const nthResultPtr nthResFun = *xpNthComb;

        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
            pool.push(std::cref(myFunParCombPerm), std::ref(parMat),
                      std::cref(v), z, n, m, step, nextStep, std::cref(freqs), myFun);
            
            SetStartZ(n, m, lower, stepSize, lowerMpz, IsRep,
                      IsComb, IsMult, IsGmp, myReps, freqs, z, nthResFun);
        }
        
        pool.push(std::cref(myFunParCombPerm), std::ref(parMat),
                  std::cref(v), z, n, m, step, nRows, std::cref(freqs), myFun);
        
        pool.join();
    } else {
        Rcpp::XPtr<combPermResPtr<typeRcpp, T>> xpFunCoPeResPtr = 
            putCombResPtrInXPtr<typeRcpp, T>(IsComb, IsMult, IsRep);
        const combPermResPtr<typeRcpp, T> myFunCombPerm = *xpFunCoPeResPtr;
        myFunCombPerm(matRcpp, v, z, n, m, 0, nRows, freqs, myFun);
    }
}

bool CheckConstrnd( SEXP f1, SEXP f2, SEXP Rtarget) {
    // No need to check myType, as we have already done 
    // so in CheckStdRet. Same goes for IsFactor. 
    bool result = !Rf_isNull(f1) && !Rf_isNull(f2) && !Rf_isNull(Rtarget);
    
    if (result) {
        if (!Rf_isString(f1))
            Rcpp::stop("constraintFun must be passed as a character");
        
        if (!Rf_isString(f2))
            Rcpp::stop("comparisonFun must be passed as a character");
    }
    
    return result;
}

// [[Rcpp::export]]
bool CheckStdRet(SEXP Rv, SEXP f1, SEXP f2, 
                 SEXP Rtarget, bool IsFactor, SEXP RKeepRes) {
    
    if (Rf_isNull(f1)) {
        return true;
    } else {
        if (IsFactor) {
            return true;
        } else {
            VecType myType = VecType::Integer;
            SetClass(myType, Rv);
            
            if (myType > VecType::Logical) {
                return true;
            } else {
                if (!Rf_isNull(f2) && !Rf_isNull(Rtarget)) {
                    return false;  // This is a constrained output
                } else if (Rf_isNull(f2) && Rf_isNull(Rtarget)) {
                    // This is applying a constrained func only
                    if (Rf_isNull(RKeepRes)) {
                        return false;
                    } else {
                        bool keepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");
                        
                        if (keepRes) {
                            return false;
                        } else {
                            return true;
                        }
                    }
                } else {
                    // This means the input is non-sensible... std return only
                    return true;
                }
            }
        }
    }
}

// [[Rcpp::export]]
SEXP CombinatoricsCnstrt(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow,
                         SEXP Rhigh, SEXP f1, SEXP f2, SEXP Rtarget, bool IsComb, 
                         SEXP RKeepRes, SEXP Rparallel, SEXP RnThreads,
                         int maxThreads, SEXP Rtolerance) {
    
    int n, m = 0, lenFreqs = 0, nRows = 0;
    bool IsMult = false;
    VecType myType = VecType::Integer;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqs;
    
    bool KeepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    
    SetClass(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    
    const bool IsConstrained = CheckConstrnd(f1, f2, Rtarget);

    if (IsConstrained) {
        for (int i = (vNum.size() - 1); i >= 0; --i)
            if (Rcpp::NumericVector::is_na(vNum[i]))
                vNum.erase(vNum.begin() + i);

        n = vNum.size();
    }
    
    SetFreqsAndM(RFreqs, IsMult, myReps, IsRep, lenFreqs, freqs, Rm, n, m);
    const std::string mainFun = Rcpp::as<std::string>(f1);
    
    if (mainFun != "prod" && mainFun != "sum" && mainFun != "mean"
            && mainFun != "max" && mainFun != "min") {
        Rcpp::stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
    }
    
    bool IsBetweenComp = false;
    std::vector<std::string> compFunVec;
    std::vector<double> targetVals;

    bool IncludeZero = false;
    PartitionType PartType = PartitionType::NotPartition;
    distinctType distinctTest;

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which cause undefined behavior
    std::vector<int> targetIntVals;
    double tolerance = 0;

    Rcpp::XPtr<funcPtr<double>> xpFunDbl = putFunPtrInXPtr<double>(mainFun);
    const funcPtr<double> funDbl = *xpFunDbl;
    
    if (IsConstrained) {                // numOnly = true, checkWhole = false, negPoss = true
        CleanConvert::convertVector(Rtarget, targetVals, "limitConstraints", true, false, true);
        compFunVec = Rcpp::as<std::vector<std::string>>(f2);
        ConstraintSetup(compFunVec, targetVals, IsBetweenComp);

        if (myType == VecType::Integer)
            if (!CheckIsInteger(mainFun, n, m, vNum, targetVals, funDbl, true))
                myType = VecType::Numeric;

        AdjustTargetVals(n, myType, targetVals, targetIntVals,
                         Rtolerance, compFunVec, tolerance, mainFun, vNum);
        
        GetPartitionCase(compFunVec, vNum, mainFun, targetVals[0], PartType,
                         distinctTest, Rlow, myReps, n, m, tolerance, IsMult, IsRep);
    }
    
    std::vector<int> startZ(m);
    double computedRows = 0;

    if (PartType > PartitionType::PartGeneral) {
        // vNum and myReps were sorted in GetPartitionCase
        IncludeZero = (vNum.front() == 0);
        int targetInt = static_cast<int>(targetVals[0]);

        SetStartPartitionZ(PartType, distinctTest, startZ,
                           myReps, targetInt, n, m, IncludeZero);
        computedRows = GetComputedPartsComps(startZ, PartType, targetInt, m,
                                             IsComb, IncludeZero, Rf_isNull(Rm));
    } else {
        computedRows = GetComputedRows(IsMult, IsComb, IsRep, n,
                                       m, Rm, lenFreqs, freqs, myReps);
    }
    
    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }

    double lower = 0, upper = 0;
    bool bLower = false, bUpper = false;
    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);

    mpz_init(lowerMpz[0]); mpz_init(upperMpz[0]);
    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz.get(), upperMpz.get(), computedRowsMpz, computedRows);

    Rcpp::XPtr<nthResultPtr> xpComb = putNthResPtrInXPtr(IsComb, IsMult, IsRep, IsGmp);
    const nthResultPtr nthResFun = *xpComb;

    if (PartType <= PartitionType::PartGeneral) {
        SetStartZ(n, m, lower, 0, lowerMpz[0], IsRep, IsComb,
                  IsMult, IsGmp, myReps, freqs, startZ, nthResFun);
    }
    
    double userNumRows = 0;
    const bool IsGenCnstrd = (IsConstrained && PartType <= PartitionType::PartGeneral);
    SetNumResults(IsGmp, bLower, bUpper, IsGenCnstrd, upperMpz.get(), lowerMpz.get(),
                  lower, upper, computedRows, computedRowsMpz, nRows, userNumRows);

    int nThreads = 1;
    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows, myType, nThreads, RnThreads, limit);

    if (IsConstrained) {
        if (myType == VecType::Integer) {
            return ConstraintReturn<Rcpp::IntegerMatrix>(n, m, mainFun, compFunVec, targetIntVals, tolerance,
                                                         vInt, bLower, lower, IsGmp, bUpper, upper, computedRows,
                                                         lowerMpz[0], upperMpz[0], computedRowsMpz, IsRep, nRows,
                                                         KeepRes, startZ, IsMult, IsComb, Rm, userNumRows, myReps,
                                                         freqs, IsBetweenComp, PartType, distinctTest.getAll);
        } else {
            return ConstraintReturn<Rcpp::NumericMatrix>(n, m, mainFun, compFunVec, targetVals, tolerance,
                                                         vNum, bLower, lower, IsGmp, bUpper, upper, computedRows,
                                                         lowerMpz[0], upperMpz[0], computedRowsMpz, IsRep, nRows,
                                                         KeepRes, startZ, IsMult, IsComb, Rm, userNumRows, myReps,
                                                         freqs, IsBetweenComp, PartType, distinctTest.getAll);
        }
    } else {
        int nCols = m + 1;

        if (myType == VecType::Integer)
            if (!CheckIsInteger(mainFun, n, m, vNum, vNum, funDbl))
                myType = VecType::Numeric;

        if (myType == VecType::Integer) {
            Rcpp::XPtr<funcPtr<int>> xpFunInt = putFunPtrInXPtr<int>(mainFun);
            const funcPtr<int> funInt = *xpFunInt;

            Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(nRows, nCols);
            MasterResRet(matInt, vInt, funInt, n, m, IsRep, IsComb, IsMult, IsGmp, freqs,
                         startZ, myReps, lower, lowerMpz[0], nRows, nThreads, Parallel);
            return matInt;
        } else {
            Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(nRows, nCols);
            MasterResRet(matNum, vNum, funDbl, n, m, IsRep, IsComb, IsMult, IsGmp, freqs,
                         startZ, myReps, lower, lowerMpz[0], nRows, nThreads, Parallel);
            return matNum;
        }
    }
}
