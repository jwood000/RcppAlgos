#include "CombinationResults.h"
#include "PermutationResults.h"
#include "PartitionsMaster.h"
#include "ConstraintsGeneral.h"
#include "GmpCombPermUtils.h"
#include "PartitionsCounts.h"
#include "RMatrix.h"
#include <RcppThread.h>

template <typename T, typename U>
using combPermResPtr = void (*const)(T &matRcpp, const std::vector<U> &v,
                             std::vector<int> z, int n, int m, int strt, int nRows,
                             const std::vector<int> &freqs, funcPtr<U> myFun);

template <typename T, typename U>
Rcpp::XPtr<combPermResPtr<T, U>> putCombResPtrInXPtr(bool IsComb, bool IsMult, bool IsRep) {
    if (IsComb) {
        if (IsMult)
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&MultisetComboResult)));
        else if (IsRep)
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&ComboGenResRep)));
        else
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&ComboGenResNoRep)));
    } else {
        if (IsMult)
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&MultisetPermRes)));
        else if (IsRep)
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&PermuteGenResRep)));
        else
            return(Rcpp::XPtr<combPermResPtr<T, U>>(new combPermResPtr<T, U>(&PermuteGenResNoRep)));
    }
};

// This is called when we can't easily produce a (loose) monotonic sequence overall,
// and we must generate and test every possible combination/permutation. This occurs
// when we are using "prod" and we have negative numbers involved. We also call this
// when lower is invoked implying that we are testing a specific range.
template <typename typeRcpp, typename T>
typeRcpp SpecCaseRet(int n, int m, std::vector<T> v, bool IsRep, int nRows, 
                     bool keepRes, std::vector<int> z, double lower, std::string mainFun, 
                     bool IsMult, double computedRows, std::vector<std::string> compFunVec,
                     std::vector<T> targetVals, bool IsComb, std::vector<int> myReps,
                     std::vector<int> freqs, bool bLower, double userRows) {
    
    if (!bLower) {
        if (computedRows > std::numeric_limits<int>::max())
            Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
        
        nRows = static_cast<int>(computedRows);
    }
    
    std::vector<T> rowVec(m);
    std::vector<int> indexMatch;
    indexMatch.reserve(nRows);
    
    Rcpp::XPtr<funcPtr<T>> xpFun = putFunPtrInXPtr<T>(mainFun);
    funcPtr<T> myFun = *xpFun;
    typeRcpp matRes = Rcpp::no_init_matrix(nRows, m + 1);
    
    // We pass keepRes = true (second true) as we need the results to determine which
    // results are within the constraint. The actual value of keepRes is utilized
    // below for the return matrix. The variable permNonTrivial, has no affect when
    // keepRes = true, so we pass it arbitrarily as true.
    Rcpp::XPtr<combPermResPtr<typeRcpp, T>> xpFunCoPeResPtr = 
        putCombResPtrInXPtr<typeRcpp, T>(IsComb, IsMult, IsRep);
    const combPermResPtr<typeRcpp, T> myFunCombPerm = *xpFunCoPeResPtr;
    myFunCombPerm(matRes, v, z, n, m, 0, nRows, freqs, myFun);
    
    Rcpp::XPtr<compPtr<T>> xpComp = putCompPtrInXPtr<T>(compFunVec[0]);
    compPtr<T> myComp = *xpComp;
    
    Rcpp::XPtr<compPtr<T>> xpComp2 = xpComp;
    compPtr<T> myComp2;
    
    if (compFunVec.size() == 1) {
        for (int i = 0; i < nRows; ++i) {
            const T testVal = matRes(i, m);
            
            if (myComp(testVal, targetVals))
                indexMatch.push_back(i);
        }
    } else {
        xpComp2 = putCompPtrInXPtr<T>(compFunVec[1]);
        myComp2 = *xpComp2;
        std::vector<T> targetVals2 = targetVals;
        targetVals2.erase(targetVals2.begin());
        
        for (int i = 0; i < nRows; ++i) {
            const T testVal = matRes(i, m);
            
            if (myComp(testVal, targetVals) || myComp2(testVal, targetVals2))
                indexMatch.push_back(i);
        }
    }
    
    const int numCols = (keepRes) ? (m + 1) : m;
    const int numMatches = indexMatch.size();
    
    if (bLower)
        nRows = numMatches;
    else
        nRows  = (numMatches > userRows) ? userRows : numMatches;
    
    typeRcpp returnMatrix = Rcpp::no_init_matrix(nRows, numCols);
    const int lastCol = keepRes ? (m + 1) : m;
    
    for (int i = 0; i < nRows; ++i)
        for (int j = 0; j < lastCol; ++j)
            returnMatrix(i, j) = matRes(indexMatch[i], j);
    
    return returnMatrix;
}

template <typename typeRcpp, typename typeVector>
typeRcpp ConstraintReturn(int n, int m, const std::string &mainFun, const std::vector<std::string> &compFunVec,
                          const std::vector<typeVector> &targetVals, double tolerance, std::vector<typeVector> &v,
                          bool bLower, double lower, bool IsGmp, bool bUpper, double upper, double computedRows,
                          mpz_t lowerMpz, mpz_t upperMpz, mpz_t computedRowsMpz, bool IsRep, int nRows, 
                          bool KeepRes, std::vector<int> &z, bool IsMult, bool IsComb, const SEXP &Rm, 
                          double userNumRows, std::vector<int> &myReps, const std::vector<int> &freqs, 
                          bool IsBetweenComp, PartitionType PartType, bool distGetAll) {
    
    const bool SpecialCase = CheckSpecialCase(n, bLower, mainFun, v);
    
    if (SpecialCase) {
        return SpecCaseRet<typeRcpp>(n, m, v, IsRep, nRows, KeepRes, z, lower, mainFun, 
                                     IsMult, computedRows, compFunVec, targetVals, IsComb,
                                     myReps, freqs, bLower, userNumRows);
    }
    
    bool bUserRows = bUpper;
    
    if (PartType > PartitionType::NotPartition) {
        std::vector<int64_t> v64(v.cbegin(), v.cend());
        int64_t target64 = static_cast<int64_t>(targetVals[0]);
        
        return Partitions::GeneralPartitions<typeRcpp>(v64, z, myReps, PartType, target64, n, m, IsRep,
                                                       IsMult, userNumRows, IsComb, KeepRes, bUserRows,
                                                       Rf_isNull(Rm), distGetAll);
    }
    
    return CombinatoricsConstraints<typeRcpp>(n, m, v, IsRep, mainFun, compFunVec, targetVals, userNumRows,
                                              IsComb, KeepRes, myReps, IsMult, bUserRows, IsBetweenComp);
}

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
            SetType(myType, Rv);
            
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
    
    SetType(myType, Rv);
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
