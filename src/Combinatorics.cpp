#include "ConstraintsMaster.h"
#include "CleanConvert.h"
#include "NthResult.h"
#include "GmpCombPermUtils.h"
#include "RMatrix.h"
#include <RcppThread.h>

// in R console: print(sqrt(.Machine$double.eps), digits = 16)
// [1] 0.00000001490116119384766
// Which is also 2^(-26)
constexpr double defaultTolerance = 0.00000001490116119384766;

// [[Rcpp::export]]
int cpp11GetNumThreads() {
    return std::thread::hardware_concurrency();
}

// N.B. Passing signed int to function expecting std::size_t is well defined
template <typename typeRcpp, typename typeVector>
void SerialReturn(int n, int r, typeVector v, bool IsRep, int nRows, bool IsComb,
                  const std::vector<int> &myReps, const std::vector<int> &freqs, 
                  std::vector<int> z, bool permNonTriv, bool IsMultiset, bool keepRes,
                  typeRcpp &matRcpp, int count, std::size_t phaseOne) {
    if (IsComb) {
        if (IsMultiset)
            MultisetCombination(n, r, v, myReps, freqs, count, nRows, z, matRcpp);
        else
            ComboGeneral(n, r, v, IsRep, count, nRows, z, matRcpp);
    } else {
        if (IsMultiset) 
            MultisetPermutation(n, r, v, nRows, z, count, matRcpp);
        else if (permNonTriv)
            PermuteGeneral(n, r, v, IsRep, nRows, z, count, matRcpp);
        else
            PermuteSerialDriver(n, r, v, IsRep, nRows, phaseOne, z, matRcpp);
    }
}

template <typename typeVector>
void ApplyFunction(int n, int r, typeVector sexpVec, bool IsRep, int nRows, bool IsComb,
                   std::vector<int> myReps, Rcpp::List &myList, std::vector<int> freqs,
                   std::vector<int> z, bool IsMultiset, SEXP sexpFun, SEXP rho, int count) {
    if (IsComb) {
        if (IsMultiset)
            MultisetComboApplyFun(n, r, sexpVec, myReps, freqs, nRows, z, count, sexpFun, rho, myList);
        else
            ComboGeneralApplyFun(n, r, sexpVec, IsRep, count, nRows, z, sexpFun, rho, myList);
    } else {
        PermutationApplyFun(n, r, sexpVec, IsRep, nRows, IsMultiset, z, count, sexpFun, rho, myList);
    }
}

// Check if our function operating on the rows of our matrix can possibly produce elements
// greater than std::numeric_limits<int>::max(). We need a NumericMatrix in this case. We also need to check
// if our function is the mean as this can produce non integral values.
bool CheckIsInteger(const std::string &funPass, std::size_t uM, int n,
                    std::vector<double> vNum, std::vector<double> targetVals,
                    funcPtr<double> myFunDbl, bool checkLim = false) {
    
    if (funPass == "mean")
        return false;
    
    std::vector<double> vAbs;
    
    for (int i = 0; i < n; ++i)
        vAbs.push_back(std::abs(vNum[i]));
    
    double vecMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
    const std::vector<double> rowVec(uM, vecMax);
    double testIfInt = myFunDbl(rowVec, uM);
    
    if (testIfInt > std::numeric_limits<int>::max())
        return false;
    
    if (checkLim) {
        vAbs.clear();
        
        for (std::size_t i = 0; i < targetVals.size(); ++i) {
            if (static_cast<int64_t>(targetVals[i]) != targetVals[i])
                return false;
            else
                vAbs.push_back(std::abs(targetVals[i]));
        }
        
        double vecMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
        
        if (vecMax > std::numeric_limits<int>::max())
            return false;
    }
    
    return true;
}

void GetStartZ(int n, int r, double &lower, int stepSize, mpz_t &lowerMpz, bool IsRep,
               bool IsComb, bool IsMultiset, bool isGmp, std::vector<int> &myReps,
               std::vector<int> freqsExpanded, std::vector<int> &startZ, nthResutlPtr nthResFun) {
    
    if (isGmp) {
        mpz_add_ui(lowerMpz, lowerMpz, stepSize);
    } else {
        lower += stepSize;
    }
    
    startZ = nthResFun(n, r, lower, lowerMpz, myReps);
    
    if (!IsComb) {
        if (IsMultiset) {
            for (std::size_t j = 0; j < startZ.size(); ++j) {
                for (std::size_t i = 0; i < freqsExpanded.size(); ++i) {
                    if (freqsExpanded[i] == startZ[j]) {
                        freqsExpanded.erase(freqsExpanded.begin() + i);
                        break;
                    }
                }
            }
            
            for (std::size_t i = 0; i < freqsExpanded.size(); ++i)
                startZ.push_back(freqsExpanded[i]);
        } else if (!IsRep) {
            if (r < n) {
                for (int i = 0; i < n; ++i) {
                    bool bExist = false;
                    for (std::size_t j = 0; j < startZ.size(); ++j) {
                        if (startZ[j] == i) {
                            bExist = true;
                            break;
                        }
                    }
                    if (!bExist)
                        startZ.push_back(i);
                }
            }
        }
    }
}

template <typename typeRcpp, typename typeElem>
void MasterReturn(int n, int r, std::vector<typeElem> v, bool IsRep, bool IsComb, std::vector<int> myReps,
                  std::vector<int> freqsExpanded, bool IsMult, std::vector<int> startZ, bool permNonTriv,
                  funcPtr<typeElem> myFun, bool keepRes, bool IsGmp, double lower, mpz_t &lowerMpz, int nRows,
                  typeRcpp &matRcpp, int nThreads, bool Parallel, std::size_t phaseOne, nthResutlPtr nthResFun) {
    
    bool generalReturn = IsComb || IsMult || permNonTriv || keepRes || n == 1;
    const std::size_t uRowN = nRows;
    
    if (!generalReturn && phaseOne > uRowN)
        generalReturn = true;
    
    if (Parallel) {
        RcppParallel::RMatrix<typeElem> parMat(matRcpp);
        
        if (generalReturn) {
            RcppThread::ThreadPool pool(nThreads);
            const int stepSize = nRows / nThreads;
            int nextStep = stepSize;
            int step = 0;
            
            for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<typeElem>, typeElem>),
                          n, r, v, IsRep, nextStep, IsComb, myReps, freqsExpanded, startZ,
                          generalReturn, IsMult, myFun, keepRes, std::ref(parMat), step, phaseOne);
                
                GetStartZ(n, r, lower, stepSize, lowerMpz, IsRep, 
                          IsComb, IsMult, IsGmp, myReps, freqsExpanded, startZ, nthResFun);
            }
            
            pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<typeElem>, typeElem>),
                      n, r, v, IsRep, nRows, IsComb, myReps, freqsExpanded, startZ,
                      generalReturn, IsMult, myFun, keepRes, std::ref(parMat), step, phaseOne);

            pool.join();
        } else {
            PermuteParallel(n, r, v, IsRep, uRowN, phaseOne, startZ, parMat, nThreads);
        }
    } else {
        GeneralReturn(n, r, v, IsRep, nRows, IsComb, myReps, freqsExpanded, 
                      startZ, generalReturn, IsMult, myFun, keepRes, matRcpp, 0, phaseOne);
    }
}

// [[Rcpp::export]]
SEXP CombinatoricsRcpp(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow,
                       SEXP Rhigh, SEXP f1, SEXP f2, SEXP Rtarget, bool IsComb, 
                       SEXP RKeepRes, bool IsFactor, bool IsCount, SEXP stdFun, 
                       SEXP myEnv, SEXP Rparallel, SEXP RNumThreads, int maxThreads,
                       SEXP Rtolerance) {
    
    int n, m = 0, lenFreqs = 0, nRows = 0;
    bool IsLogical, IsCharacter, IsMultiset, IsInteger, IsComplex, IsRaw;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqsExpanded;
    Rcpp::CharacterVector rcppChar;
    Rcpp::ComplexVector rcppCplx;
    Rcpp::RawVector rcppRaw;
    
    bool keepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRepetition = CleanConvert::convertLogical(RisRep, "repetition");
    SetClass(IsCharacter, IsLogical, IsInteger, IsComplex, IsRaw, Rv);
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        myReps.push_back(1);
    } else {
        IsRepetition = false;
        CleanConvert::convertVector(RFreqs, myReps, "freqs");
        lenFreqs = static_cast<int>(myReps.size());
        bool allOne = std::all_of(myReps.cbegin(), myReps.cend(), 
                                    [](int v_i) {return v_i == 1;});
        if (allOne) {
            IsMultiset = false;
            freqsExpanded = myReps;
        } else {
            IsMultiset = true;
            
            for (int i = 0; i < lenFreqs; ++i)
                for (int j = 0; j < myReps[i]; ++j)
                    freqsExpanded.push_back(i);
        }
    }
    
    const bool mIsNull = Rf_isNull(Rm);
    
    if (mIsNull) {
        if (!freqsExpanded.empty())
            m = freqsExpanded.size();
    } else {
        if (Rf_length(Rm) > 1)
            Rcpp::stop("length of m must be 1");
        
        CleanConvert::convertPrimitive(Rm, m, "m");
    }
    
    SetValues(IsCharacter, IsLogical, IsInteger, IsComplex, 
              IsRaw, rcppChar, vInt, vNum, rcppCplx, rcppRaw, n, Rv);
    
    bool IsConstrained = false; 

    if (IsFactor)
        keepRes = IsComplex = IsCharacter = IsInteger = false;

    if (!(IsLogical || IsCharacter || IsComplex || IsRaw
              || Rf_isNull(f1) || Rf_isNull(f2) || Rf_isNull(Rtarget))) {
        if (!Rf_isString(f1))
            Rcpp::stop("constraintFun must be passed as a character");

        if (!Rf_isString(f2))
            Rcpp::stop("comparisonFun must be passed as a character");
        
        IsConstrained = true;
    }

    if (IsConstrained) {
        for (int i = (vNum.size() - 1); i >= 0; --i)
            if (Rcpp::NumericVector::is_na(vNum[i]))
                vNum.erase(vNum.begin() + i);

        n = vNum.size();
    }

    if (!IsMultiset && mIsNull && freqsExpanded.empty()) {m = n;}
    const double computedRows = GetComputedRows(IsMultiset, IsComb, IsRepetition, n,
                                                m, Rm, lenFreqs, freqsExpanded, myReps);

    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);

    if (IsGmp) {
        GetComputedRowMpz(computedRowMpz, IsMultiset, IsComb,
                          IsRepetition, n, m, Rm, freqsExpanded, myReps);
    }

    double lower = 0, upper = 0;
    bool bLower = false, bUpper = false;
    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);

    mpz_init(lowerMpz[0]); mpz_init(upperMpz[0]);
    SetBounds(IsCount, Rlow, Rhigh, IsGmp, bLower, bUpper,
              lower, upper, lowerMpz.get(), upperMpz.get());

    // N.B. for the lower analogs we are using >= as we have subtracted one
    // Also, we are postponing this lower is not given and IsConstrained is true
    // as there is a possibility that are computed number of results is not
    // correct. This can happen in cases where we are finding partitions.
    if (!(IsConstrained && !bLower))
        CheckBounds(IsGmp, lower, upper, computedRows, lowerMpz[0], upperMpz[0], computedRowMpz);

    if (IsCount)
        return GetCount(IsGmp, computedRowMpz, computedRows);

    std::vector<int> startZ(m);
    bool permNonTriv = false;
    double dblLower = lower;
    if (!IsGmp) mpz_set_d(lowerMpz[0], dblLower);

    Rcpp::XPtr<nthResutlPtr> xpComb = putNthResPtrInXPtr(IsMultiset, IsRepetition, IsGmp, IsComb);
    nthResutlPtr nthResFun = *xpComb;

    if (bLower && mpz_cmp_ui(lowerMpz[0], 0) > 0) {
        GetStartZ(n, m, lower, 0, lowerMpz[0], IsRepetition, IsComb,
                  IsMultiset, IsGmp, myReps, freqsExpanded, startZ, nthResFun);
        if (!IsComb) permNonTriv = true;
    } else {
        if (IsComb) {
            if (IsMultiset)
                startZ.assign(freqsExpanded.cbegin(), freqsExpanded.cbegin() + m);
            else if (IsRepetition)
                std::fill(startZ.begin(), startZ.end(), 0);
            else
                std::iota(startZ.begin(), startZ.end(), 0);
        } else {
            if (IsMultiset) {
                startZ = freqsExpanded;
            } else if (IsRepetition) {
                std::fill(startZ.begin(), startZ.end(), 0);
            } else {
                startZ.resize(n);
                std::iota(startZ.begin(), startZ.end(), 0);
            }
        }
    }

    double userNumRows = 0;
    SetNumResults(IsGmp, bLower, bUpper, IsConstrained, permNonTriv, upperMpz.get(),
                  lowerMpz.get(), lower, upper, computedRows, computedRowMpz, nRows, userNumRows);

    const std::size_t uM = m;
    int nThreads = 1;
    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows, IsCharacter, nThreads, RNumThreads, limit);

    if (IsConstrained) {
        std::vector<double> targetVals;      // numOnly = true, checkWhole = false, negPoss = true
        CleanConvert::convertVector(Rtarget, targetVals, "limitConstraints", true, false, true);

        if (targetVals.size() > 2)
            Rcpp::stop("there cannot be more than 2 limitConstraints");
        else if (targetVals.size() == 2 && targetVals[0] == targetVals[1])
            Rcpp::stop("The limitConstraints must be different");

        const std::string mainFun = Rcpp::as<std::string>(f1);
        if (mainFun != "prod" && mainFun != "sum" && mainFun != "mean"
                && mainFun != "max" && mainFun != "min") {
            Rcpp::stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
        }

        std::vector<std::string>::const_iterator itComp;
        std::vector<std::string> compFunVec = Rcpp::as<std::vector<std::string>>(f2);

        if (compFunVec.size() > 2)
            Rcpp::stop("there cannot be more than 2 comparison operators");

        for (std::size_t i = 0; i < compFunVec.size(); ++i) {
            itComp = std::find(compForms.cbegin(), compForms.cend(), compFunVec[i]);

            if (itComp == compForms.end()) {
                Rcpp::stop("comparison operators must be one of the following: "
                               "'>', '>=', '<', '<=', or '=='");
            }

            int myIndex = std::distance(compForms.cbegin(), itComp);

            // The first 5 are "standard" whereas the 6th and 7th
            // are written with the equality first. Converting
            // them here makes it easier to deal with later.
            if (myIndex > 4) {myIndex -= 3;}
            compFunVec[i] = compForms[myIndex];
        }

        bool between = false;

        if (compFunVec.size() == 2) {
            if (targetVals.size() == 1) {
                compFunVec.pop_back();
            } else {
                if (compFunVec[0] == "==" || compFunVec[1] == "==") {
                    Rcpp::stop("If comparing against two limitConstraints, the "
                                "equality comparisonFun (i.e. '==') cannot be used. "
                                "Instead, use '>=' or '<='.");
                }

                if (compFunVec[0].substr(0, 1) == compFunVec[1].substr(0, 1)) {
                    Rcpp::stop("Cannot have two 'less than' comparisonFuns or two 'greater than' "
                          "comparisonFuns (E.g. c('<', '<=') is not allowed).");
                }

                // The two cases below are for when we are looking for all combs/perms such that when
                // myFun is applied, the result is between 2 values. These comparisons are defined in
                // compSpecial in ConstraintsMaster.h. If we look at the definitions of these comparison
                // functions in ConstraintsUtils.h, we see the following trend for all 4 cases:
                //
                // bool greaterLess(stdType x, const std::vector<stdType> &y) {return (x < y[0]) && (x > y[1]);}
                //
                // That is, we need the maixmum value in targetVals to be the first value and the second
                // value needs to be the minimum. At this point, the constraint algorithm will be
                // identical to when comparisonFun = "==" (i.e. we allow the algorithm to continue
                // while the results are less than (or equal to in cases where strict inequalities are
                // enforced) the target value and stop once the result exceeds that value).

                if (compFunVec[0].substr(0, 1) == ">" && std::min(targetVals[0], targetVals[1]) == targetVals[0]) {
                    compFunVec[0] = compFunVec[0] + "," + compFunVec[1];
                    between = true;
                } else if (compFunVec[0].substr(0, 1) == "<" && std::max(targetVals[0], targetVals[1]) == targetVals[0]) {
                    compFunVec[0] = compFunVec[1] + "," + compFunVec[0];
                    between = true;
                }

                if (between) {
                    compFunVec.pop_back();

                    if (std::max(targetVals[0], targetVals[1]) == targetVals[1])
                        std::swap(targetVals[0], targetVals[1]);
                }
            }
        } else {
            if (targetVals.size() == 2)
                targetVals.pop_back();
        }

        Rcpp::XPtr<funcPtr<double>> xpFunDbl = putFunPtrInXPtr<double>(mainFun);
        funcPtr<double> myFunDbl = *xpFunDbl;

        IsInteger = (IsInteger) && CheckIsInteger(mainFun, uM, n, vNum, targetVals, myFunDbl, true);
        // Must be defined inside IsInteger check as targetVals could be outside integer data type range
        std::vector<int> targetIntVals;
        double tolerance = 0;

        if (IsInteger) {
            targetIntVals.assign(targetVals.cbegin(), targetVals.cend());
        } else {
            // We first check if we are getting double precision. If so, for the
            // non-strict inequalities, we must alter the limit by some epsilon:
            //
            //                 x <= y   --->>>   x <= y + e
            //                 x >= y   --->>>   x >= y - e
            //
            // Equality is a bit tricky as we need to check whether the absolute
            // value of the difference is less than epsilon. As a result, we
            // can't alter the limit with one alteration. Observe:
            //
            //          x == y  --->>>  |x - y| <= e , which gives:
            //
            //                      -e <= x - y <= e
            //
            //                    1.     x >= y - e
            //                    2.     x <= y + e
            //
            // As a result, we must define a specialized equality check for double
            // precision. It is 'equalDbl' and can be found in ConstraintsUtils.h

            if (!IsInteger) {
                if (Rf_isNull(Rtolerance)) {
                    bool IsWhole = true;

                    for (int i = 0; i < n && IsWhole; ++i)
                        if (static_cast<int64_t>(vNum[i]) != vNum[i])
                            IsWhole = false;

                    for (std::size_t i = 0; i < targetVals.size() && IsWhole; ++i)
                        if (static_cast<int64_t>(targetVals[i]) != targetVals[i])
                            IsWhole = false;

                    tolerance = (IsWhole && mainFun != "mean") ? 0 : defaultTolerance;
                } else {
                    CleanConvert::convertPrimitive(Rtolerance, tolerance, "tolerance", true, false, false, true);
                }

                auto itComp = std::find(compSpecial.cbegin(),
                                        compSpecial.cend(), compFunVec[0]);

                if (compFunVec[0] == "==") {
                    targetVals.push_back(targetVals[0] - tolerance);
                    targetVals[0] += tolerance;
                } else if (itComp != compSpecial.end()) {
                    targetVals[0] += tolerance;
                    targetVals[1] -= tolerance;
                } else if (compFunVec[0] == "<=") {
                    targetVals[0] += tolerance;
                } else if (compFunVec[0] == ">=") {
                    targetVals[0] -= tolerance;
                }

                if (compFunVec.size() > 1) {
                    if (compFunVec[1] == "<=") {
                        targetVals[1] += tolerance;
                    } else if (compFunVec[1] == ">=") {
                        targetVals[1] -= tolerance;
                    }
                }
            }
        }

        bool SpecialCase = false;

        // If bLower, the user is looking to test a particular range. Otherwise, the constraint algo
        // will simply return (upper - lower) # of combinations/permutations that meet the criteria
        if (bLower) {
            SpecialCase = true;
        } else if (mainFun == "prod") {
            for (int i = 0; i < n; ++i) {
                if (vNum[i] < 0) {
                    SpecialCase = true;
                    break;
                }
            }
        }

        if (SpecialCase) {
            CheckBounds(IsGmp, lower, upper, computedRows, lowerMpz[0], upperMpz[0], computedRowMpz);

            if (IsInteger) {
                return SpecCaseRet<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, nRows, keepRes, startZ, lower,
                                                        mainFun, IsMultiset, computedRows, compFunVec, targetIntVals, IsComb,
                                                        myReps, freqsExpanded, bLower, userNumRows);
            } else {
                return SpecCaseRet<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes, startZ, lower,
                                                        mainFun, IsMultiset, computedRows, compFunVec, targetVals, IsComb,
                                                        myReps, freqsExpanded, bLower, userNumRows);
            }
        }

        bool bUserRows = bUpper;

        if (compFunVec[0] == "==" && mainFun == "sum" && n > 1 && m > 1) {
            std::vector<double> pTest(vNum.cbegin(), vNum.cend());
            std::sort(pTest.begin(), pTest.end());
            const double tarDiff = pTest[1] - pTest[0];

            if (static_cast<int64_t>(pTest[0]) == pTest[0]) {
                bool PartitionCase = true;

                for (int i = 1; i < n; ++i) {
                    const double testDiff = pTest[i] - pTest[i - 1];

                    if (std::abs(testDiff - tarDiff) > tolerance
                            || static_cast<int64_t>(pTest[i]) != pTest[i]) {
                        PartitionCase = false;
                        break;
                    }
                }

                if (PartitionCase) {
                    int64_t target64 = static_cast<int64_t>(targetVals[0]);
                    std::vector<int64_t> v64(vNum.cbegin(), vNum.cend());

                    if (target64 == targetVals[0]) {
                        if (IsInteger) {
                            return Partitions::GeneralPartitions<Rcpp::IntegerMatrix>(n, m, v64, target64, IsRepetition, IsMultiset,
                                                                                      myReps, userNumRows, IsComb, keepRes, bUserRows, mIsNull);
                        } else {
                            return Partitions::GeneralPartitions<Rcpp::NumericMatrix>(n, m, v64, target64, IsRepetition, IsMultiset,
                                                                                      myReps, userNumRows, IsComb, keepRes, bUserRows, mIsNull);
                        }
                    }
                }
            }
        }

        CheckBounds(IsGmp, lower, upper, computedRows, lowerMpz[0], upperMpz[0], computedRowMpz);

        if (IsInteger) {
            return CombinatoricsConstraints<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, mainFun, compFunVec, targetIntVals,
                                                                 userNumRows, IsComb, keepRes, myReps, IsMultiset, bUserRows, between);
        }

        return CombinatoricsConstraints<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, mainFun, compFunVec, targetVals,
                                                             userNumRows, IsComb, keepRes, myReps, IsMultiset, bUserRows, between);
    } else {
        bool applyFun = !Rf_isNull(stdFun) && !IsFactor;

        if (applyFun) {
            if (!Rf_isFunction(stdFun))
                Rcpp::stop("FUN must be a function!");

            Rcpp::List myList(nRows);
            SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));

            if (IsCharacter) {
                ApplyFunction(n, m, rcppChar, IsRepetition, nRows, IsComb, myReps,
                              myList, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            } else if (IsComplex) {
                ApplyFunction(n, m, rcppCplx, IsRepetition, nRows, IsComb, myReps,
                              myList, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            } else if (IsRaw) {
                ApplyFunction(n, m, rcppRaw, IsRepetition, nRows, IsComb, myReps,
                              myList, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            } else if (IsLogical || IsInteger) {
                Rcpp::IntegerVector rcppVInt(vInt.cbegin(), vInt.cend());
                ApplyFunction(n, m, rcppVInt, IsRepetition, nRows, IsComb, myReps,
                              myList, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            } else {
                Rcpp::NumericVector rcppVNum(vNum.cbegin(), vNum.cend());
                ApplyFunction(n, m, rcppVNum, IsRepetition, nRows, IsComb, myReps,
                              myList, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            }

            UNPROTECT(1);
            return myList;
        }

        // It is assumed that if user has constraintFun with no comparison
        // or limitConstraints and they have not explicitly set keepRes to
        // FALSE, then they simply want the constraintFun applied
        if (Rf_isNull(RKeepRes)) {
            if (Rf_isNull(f2) && Rf_isNull(Rtarget) && !Rf_isNull(f1))
                keepRes = !IsLogical && !IsCharacter && !IsFactor;
        } else {
            keepRes = keepRes && !Rf_isNull(f1);
        }
        
        funcPtr<double> myFunDbl = *putFunPtrInXPtr<double>("sum");
        funcPtr<int> myFunInt = *putFunPtrInXPtr<int>("sum");
        int nCol = m;

        if (keepRes) {
            const std::string mainFun = Rcpp::as<std::string>(f1);
            if (mainFun != "prod" && mainFun != "sum" && mainFun != "mean"
                    && mainFun != "max" && mainFun != "min") {
                Rcpp::stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
            }
            
            Rcpp::XPtr<funcPtr<double>> xpFunDbl = putFunPtrInXPtr<double>(mainFun);
            Rcpp::XPtr<funcPtr<int>> xpFunInt = putFunPtrInXPtr<int>(mainFun);
            myFunDbl = *xpFunDbl;
            myFunInt = *xpFunInt;

            IsInteger = (IsInteger) && CheckIsInteger(mainFun, uM, n, vNum, vNum, myFunDbl);
            ++nCol;
        }

        const std::size_t phaseOne = (!permNonTriv && !IsComb)
                                     ? ((IsRepetition) ? std::pow(static_cast<double>(n),
                                                                  static_cast<double>(m - 1))
                                                       : NumPermsNoRep(n - 1, m - 1)) : 0u;
        if (IsCharacter) {
            Rcpp::CharacterMatrix matChar = Rcpp::no_init_matrix(nRows, nCol);
            SerialReturn(n, m, rcppChar, IsRepetition, nRows, IsComb, myReps,
                         freqsExpanded, startZ, permNonTriv, IsMultiset, keepRes,
                         matChar, 0, phaseOne);
            return matChar;
        } else if (IsComplex) {
            Rcpp::ComplexMatrix matCplx(nRows, nCol);
            SerialReturn(n, m, rcppCplx, IsRepetition, nRows, IsComb, myReps,
                         freqsExpanded, startZ, permNonTriv, IsMultiset, keepRes,
                         matCplx, 0, phaseOne);
            return matCplx;
        } else if (IsRaw) {
            Rcpp::RawMatrix matRaw(nRows, nCol);
            SerialReturn(n, m, rcppRaw, IsRepetition, nRows, IsComb, myReps,
                         freqsExpanded, startZ, permNonTriv, IsMultiset, keepRes,
                         matRaw, 0, phaseOne);
            return matRaw;
        } else if (IsLogical) {
            Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(nRows, nCol);
            MasterReturn(n, m, vInt, IsRepetition, IsComb, myReps, freqsExpanded,
                         IsMultiset, startZ, permNonTriv, myFunInt, keepRes, IsGmp, lower,
                         lowerMpz[0], nRows, matBool, nThreads, Parallel, phaseOne, nthResFun);
            return matBool;
        } else if (IsFactor || IsInteger) {
            Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(nRows, nCol);
            MasterReturn(n, m, vInt, IsRepetition, IsComb, myReps, freqsExpanded,
                         IsMultiset, startZ, permNonTriv, myFunInt, keepRes, IsGmp, lower,
                         lowerMpz[0], nRows, matInt, nThreads, Parallel, phaseOne, nthResFun);
            
            if (IsFactor) {
                Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
                Rcpp::CharacterVector myClass = testFactor.attr("class");
                Rcpp::CharacterVector myLevels = testFactor.attr("levels");
                matInt.attr("class") = myClass;
                matInt.attr("levels") = myLevels;
            }

            return matInt;
        } else {
            Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(nRows, nCol);
            MasterReturn(n, m, vNum, IsRepetition, IsComb, myReps, freqsExpanded,
                         IsMultiset, startZ, permNonTriv, myFunDbl, keepRes, IsGmp, lower,
                         lowerMpz[0], nRows, matNum, nThreads, Parallel, phaseOne, nthResFun);
            return matNum;
        }
    }
}
