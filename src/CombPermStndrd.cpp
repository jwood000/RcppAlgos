#include "Combinations.h"
#include "Permutations.h"
#include "GmpCombPermUtils.h"
#include "RMatrix.h"

template <typename T, typename U>
using combPermPtr = void (*const)(T &matRcpp, const U &v, std::vector<int> z,
                          int n, int m, int strt, int nRows, const std::vector<int> &freqs);

template <typename T, typename U>
Rcpp::XPtr<combPermPtr<T, U>> putCombPtrInXPtr(bool IsComb, bool IsMult, bool IsRep, bool IsGen) {
    
    if (IsComb) {
        if (IsMult)
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&MultisetCombination)));
        else if (IsRep)
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&CombinationsRep)));
        else
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&CombinationsNoRep)));
    } else {
        if (IsMult) {
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&MultisetPermutation)));
        } else if (IsGen) {
            if (IsRep)
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteGeneralRep)));
            else
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteGeneralNoRep)));
        } else {
            if (IsRep)
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteSerialRep)));
            else
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteSerialNoRep)));
        }
    }
}

template <int RTYPE>
Rcpp::List ApplyFunction(const Rcpp::Vector<RTYPE> &v, int n, int m, bool IsRep, int nRows,
                         bool IsComb, const std::vector<int> &freqs, std::vector<int> &z, 
                         bool IsMult, SEXP stdFun, SEXP rho) {
    
    Rcpp::List myList(nRows);
    SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));
    
    if (IsComb) {
        if (IsMult)
            MultisetComboApplyFun(myList, v, z, n, m, nRows, sexpFun, rho, freqs);
        else
            ComboGeneralApplyFun(myList, v, z, n, m, IsRep, nRows, sexpFun, rho);
    } else {
        PermutationApplyFun(myList, v, z, n, m, IsRep, IsMult, nRows, sexpFun, rho);
    }
    
    UNPROTECT(1);
    return myList;
}

template <int T>
Rcpp::Matrix<T> SerialReturn(const Rcpp::Vector<T> &v, std::vector<int> &z, int n,
                             int m, int nRows, bool IsComb, bool IsRep, bool IsMult,
                             bool generalRet, const std::vector<int> &freqs) {

    Rcpp::Matrix<T> matRcpp = Rcpp::no_init_matrix(nRows, m);
    Rcpp::XPtr<combPermPtr<Rcpp::Matrix<T>, Rcpp::Vector<T>>> xpFunCoPePtr = 
        putCombPtrInXPtr<Rcpp::Matrix<T>, Rcpp::Vector<T>>(IsComb, IsMult, IsRep, generalRet);

    const combPermPtr<Rcpp::Matrix<T>, Rcpp::Vector<T>> myFunCombPerm = *xpFunCoPePtr; 
    myFunCombPerm(matRcpp, v, z, n, m, 0, nRows, freqs);
    return matRcpp;
}

template <typename typeRcpp, typename T>
void MasterReturn(typeRcpp &matRcpp, std::vector<T> v, int n, int m, bool IsRep, 
                  bool IsComb, bool IsMult, bool IsGmp, bool generalRet, const std::vector<int> &freqs,
                  std::vector<int> z, const std::vector<int> &myReps, double lower, mpz_t &lowerMpz,
                  int nRows, int nThreads, bool Parallel, int phaseOne) {

    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(matRcpp);

        if (generalRet) {
            RcppThread::ThreadPool pool(nThreads);
            const int stepSize = nRows / nThreads;
            int nextStep = stepSize;
            int step = 0;
            
            Rcpp::XPtr<combPermPtr<RcppParallel::RMatrix<T>, 
                                   std::vector<T>>> xpFunCoPePtr = 
                    putCombPtrInXPtr<RcppParallel::RMatrix<T>, 
                                     std::vector<T>>(IsComb, IsMult, IsRep, generalRet);
            
            const combPermPtr<RcppParallel::RMatrix<T>, 
                              std::vector<T>> myFunParCombPerm = *xpFunCoPePtr;
            
            Rcpp::XPtr<nthResultPtr> xpNthComb = putNthResPtrInXPtr(IsComb, IsMult, IsRep, IsGmp);
            const nthResultPtr nthResFun = *xpNthComb;

            for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                pool.push(std::cref(myFunParCombPerm), std::ref(parMat),
                          std::cref(v), z, n, m, step, nextStep, std::cref(freqs));
                
                SetStartZ(n, m, lower, stepSize, lowerMpz, IsRep,
                          IsComb, IsMult, IsGmp, myReps, freqs, z, nthResFun);
            }
            
            pool.push(std::cref(myFunParCombPerm), std::ref(parMat),
                      std::cref(v), z, n, m, step, nRows, std::cref(freqs));

            pool.join();
        } else {
            PermuteParallel(parMat, v, z, n, m, nRows, phaseOne, nThreads, IsRep);
        }
    } else {
        Rcpp::XPtr<combPermPtr<typeRcpp, std::vector<T>>> xpFunCoPePtr = putCombPtrInXPtr<typeRcpp, 
                               std::vector<T>>(IsComb, IsMult, IsRep, generalRet);
        const combPermPtr<typeRcpp, std::vector<T>> myFunCombPerm = *xpFunCoPePtr;
        myFunCombPerm(matRcpp, v, z, n, m, 0, nRows, freqs);
    }
}

// [[Rcpp::export]]
SEXP CombinatoricsCount(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, bool IsComb) {
    
    int n, m = 0, lenFreqs = 0;
    VecType myType = VecType::Integer;
    bool IsMult = false;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqs;
    
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    SetClass(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    SetFreqsAndM(RFreqs, IsMult, myReps, IsRep, lenFreqs, freqs, Rm, n, m);
    
    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep, n,
                                                m, Rm, lenFreqs, freqs, myReps);
    
    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);
    
    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult, 
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }
    
    return GetCount(IsGmp, computedRowsMpz, computedRows);
}

// [[Rcpp::export]]
SEXP CombinatoricsStndrd(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow,
                         SEXP Rhigh, bool IsComb, bool IsFactor, SEXP stdFun, 
                         SEXP myEnv, SEXP Rparallel, SEXP RNumThreads, int maxThreads) {

    int n, m = 0, lenFreqs = 0, nRows = 0;
    VecType myType = VecType::Integer;
    bool IsMult = false;

    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqs;

    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");

    SetClass(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    SetFreqsAndM(RFreqs, IsMult, myReps, IsRep, lenFreqs, freqs, Rm, n, m);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep, n,
                                                m, Rm, lenFreqs, freqs, myReps);

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

    std::vector<int> startZ(m);
    const bool permNonTriv = (!IsComb && bLower);
    
    Rcpp::XPtr<nthResultPtr> xpComb = putNthResPtrInXPtr(IsComb, IsMult, IsRep, IsGmp);
    const nthResultPtr nthResFun = *xpComb;
    
    SetStartZ(n, m, lower, 0, lowerMpz[0], IsRep, IsComb,
              IsMult, IsGmp, myReps, freqs, startZ, nthResFun);
    
    double userNumRows = 0;
    SetNumResults(IsGmp, bLower, bUpper, false, upperMpz.get(), lowerMpz.get(),
                  lower, upper, computedRows, computedRowsMpz, nRows, userNumRows);
    
    const bool applyFun = !Rf_isNull(stdFun) && !IsFactor;
    
    if (applyFun) {
        if (!Rf_isFunction(stdFun))
            Rcpp::stop("FUN must be a function!");
        
        const SEXP sexpCopy = CopyRv(Rv, vInt, vNum, myType);
        RCPP_RETURN_VECTOR(ApplyFunction, sexpCopy, n, m, IsRep, nRows,
                           IsComb, freqs, startZ, IsMult, stdFun, myEnv);
    }
    
    int nThreads = 1;
    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows, myType, nThreads, RNumThreads, limit);

    const double phaseOneDbl = (!permNonTriv && !IsComb) ?
                               ((IsRep) ? std::pow(static_cast<double>(n),
                                             static_cast<double>(m - 1)) :
                               NumPermsNoRep(n - 1, m - 1)) : 0;
    
    bool generalRet = IsComb || IsMult || permNonTriv || 
                      n == 1 || phaseOneDbl > std::numeric_limits<int>::max();
    
    if (!generalRet && phaseOneDbl > nRows) {generalRet = true;}
    const int phaseOne = (generalRet) ? 0 : static_cast<int>(phaseOneDbl);

    if (myType > VecType::Logical) {
        const SEXP sexpCopy(Rcpp::clone(Rv));
        RCPP_RETURN_VECTOR(SerialReturn, sexpCopy, startZ, n, m, 
                           nRows, IsComb, IsRep, IsMult, generalRet, freqs);
    } else if (myType == VecType::Logical) {
        Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(nRows, m);
        MasterReturn(matBool, vInt, n, m, IsRep, IsComb, IsMult,
                     IsGmp, generalRet, freqs, startZ, myReps, lower,
                     lowerMpz[0], nRows, nThreads, Parallel, phaseOne);
        return matBool;
    } else if (myType == VecType::Integer) {
        Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(nRows, m);
        MasterReturn(matInt, vInt, n, m, IsRep, IsComb, IsMult,
                     IsGmp, generalRet, freqs, startZ, myReps, lower,
                     lowerMpz[0], nRows, nThreads, Parallel, phaseOne);
        if (IsFactor) {SetFactorClass(matInt, Rv);}
        return matInt;
    } else {
        Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(nRows, m);
        MasterReturn(matNum, vNum, n, m, IsRep, IsComb, IsMult,
                     IsGmp, generalRet, freqs, startZ, myReps, lower,
                     lowerMpz[0], nRows, nThreads, Parallel, phaseOne);
        return matNum;
    }
}

