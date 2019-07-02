#include "ConstraintsMaster.h"
#include "GeneralPartitions.h"
#include "CleanConvert.h"
#include "RMatrix.h"
#include <RcppThread.h>

// [[Rcpp::export]]
int cpp11GetNumThreads() {
    return std::thread::hardware_concurrency();
}

// [[Rcpp::export]]
SEXP CombinatoricsRcpp(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow,
                       SEXP Rhigh, SEXP f1, SEXP f2, SEXP Rtarget, bool IsComb, 
                       SEXP RKeepRes, bool IsFactor, bool IsCount, SEXP stdFun, 
                       SEXP myEnv, SEXP Rparallel, SEXP RNumThreads, int maxThreads,
                       SEXP Rtolerance) {
    
    int n, m1, m2, m = 0, lenFreqs = 0, nRows = 0;
    bool IsLogical, IsCharacter, IsMultiset, IsInteger;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqsExpanded;
    Rcpp::CharacterVector rcppChar;
    
    bool keepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsRepetition = CleanConvert::convertLogical(RisRep, "repetition");
    
    switch(TYPEOF(Rv)) {
        case LGLSXP: {
            IsLogical = true;
            keepRes = IsInteger = IsCharacter = false;
            break;
        }
        case INTSXP: {
            IsInteger = true;
            IsLogical = IsCharacter = false;
            break;
        }
        case REALSXP: {
            IsLogical = IsInteger = IsCharacter = false;
            break;
        }
        case STRSXP: {
            IsCharacter = true;
            Parallel = keepRes = IsLogical = IsInteger = false;
            break;
        }
        default: {
            Rcpp::stop("Only integers, numerical, character, and factor classes are supported for v");   
        }
    }
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        myReps.push_back(1);
    } else {
        IsMultiset = true;
        IsRepetition = false;
        CleanConvert::convertVector(RFreqs, myReps, "freqs");
        lenFreqs = static_cast<int>(myReps.size());
        
        for (int i = 0; i < lenFreqs; ++i)
            for (int j = 0; j < myReps[i]; ++j)
                freqsExpanded.push_back(i);
    }
    
    if (Rf_isNull(Rm)) {
        if (IsMultiset) {
            m = freqsExpanded.size();
        } else {
            Rcpp::stop("m and freqs cannot both be NULL");
        }
    } else {
        if (Rf_length(Rm) > 1)
            Rcpp::stop("length of m must be 1");
        
        CleanConvert::convertPrimitive(Rm, m, "m");
    }
    
    if (IsCharacter) {
        rcppChar = Rcpp::as<Rcpp::CharacterVector>(Rv);
        n = rcppChar.size();
    } else if (IsLogical) {
        vInt = Rcpp::as<std::vector<int>>(Rv);
        n = vInt.size();
    } else {
        if (Rf_length(Rv) == 1) {
            int seqEnd;             // numOnly = true, checkWhole = true, negPoss = true
            CleanConvert::convertPrimitive(Rv, seqEnd, "If v is not a character and of length 1, it", true, true, true);
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            Rcpp::IntegerVector vTemp = Rcpp::seq(m1, m2);
            IsInteger = true;
            vNum = Rcpp::as<std::vector<double>>(vTemp);
        } else {
            vNum = Rcpp::as<std::vector<double>>(Rv);
        }
        
        n = vNum.size();
    }
    
    if (IsInteger) {
        for (int i = 0; i < n && IsInteger; ++i)
            if (Rcpp::NumericVector::is_na(vNum[i]))
                IsInteger = false;
        
        if (IsInteger)
            vInt.assign(vNum.cbegin(), vNum.cend());
    }
        
    bool IsConstrained = false;
    
    if (IsFactor)
        keepRes = IsConstrained = IsCharacter = IsInteger = false;
    
    if (IsLogical || IsCharacter || Rf_isNull(f1) || Rf_isNull(f2) || Rf_isNull(Rtarget)) {
        IsConstrained = false;
    } else {
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

    double computedRows = 0;
    
    if (IsMultiset) {
        if (n != lenFreqs)
            Rcpp::stop("the length of freqs must equal the length of v");
        
        if (m > static_cast<int>(freqsExpanded.size()))
            m = freqsExpanded.size();
        
        if (IsComb) {
            computedRows = MultisetCombRowNum(n, m, myReps);
        } else {
            if (Rf_isNull(Rm) || m == static_cast<int>(freqsExpanded.size()))
                computedRows = NumPermsWithRep(freqsExpanded);
            else
                computedRows = MultisetPermRowNum(n, m, myReps);
        }
    } else {
        if (IsRepetition) {
            if (IsComb)
                computedRows = NumCombsWithRep(n, m);
            else
                computedRows = std::pow(static_cast<double>(n), static_cast<double>(m));
        } else {
            if (m > n)
                Rcpp::stop("m must be less than or equal to the length of v");
            
            if (IsComb)
                computedRows = nChooseK(n, m);
            else
                computedRows = NumPermsNoRep(n, m);
        }
    }
    
    bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);
    
    if (IsGmp) {
        if (IsMultiset) {
            if (IsComb) {
                MultisetCombRowNumGmp(computedRowMpz, n, m, myReps);
            } else {
                if (Rf_isNull(Rm) || m == static_cast<int>(freqsExpanded.size()))
                    NumPermsWithRepGmp(computedRowMpz, freqsExpanded);
                else
                    MultisetPermRowNumGmp(computedRowMpz, n, m, myReps);
            }
        } else {
            if (IsRepetition) {
                if (IsComb)
                    NumCombsWithRepGmp(computedRowMpz, n, m);
                else
                    mpz_ui_pow_ui(computedRowMpz, n, m);
            } else {
                if (IsComb)
                    nChooseKGmp(computedRowMpz, n, m);
                else
                    NumPermsNoRepGmp(computedRowMpz, n, m);
            }
        }
    }

    double lower = 0, upper = 0;
    bool bLower = false, bUpper = false;
    mpz_t lowerMpz[1], upperMpz[1];
    mpz_init(lowerMpz[0]); mpz_init(upperMpz[0]);
    mpz_set_ui(lowerMpz[0], 0); mpz_set_ui(upperMpz[0], 0);
    
    if (!IsCount) {
        if (!Rf_isNull(Rlow)) {
            bLower = true;
            
            if (IsGmp) {
                createMPZArray(Rlow, lowerMpz, 1, "lower");
                mpz_sub_ui(lowerMpz[0], lowerMpz[0], 1);
            } else {                                    // numOnly = false
                CleanConvert::convertPrimitive(Rlow, lower, "lower", false);
                --lower;
            }
        }
        
        if (!Rf_isNull(Rhigh)) {
            bUpper = true;
            
            if (IsGmp) {
                createMPZArray(Rhigh, upperMpz, 1, "upper");
            } else {                                     // numOnly = false
                CleanConvert::convertPrimitive(Rhigh, upper, "upper", false);
            }
        }
    }
    
    // N.B. for the lower analogs we are using >= as we have subtracted one
    if (IsGmp) {
        if (mpz_cmp(lowerMpz[0], computedRowMpz) >= 0 || mpz_cmp(upperMpz[0], computedRowMpz) > 0)
            Rcpp::stop("bounds cannot exceed the maximum number of possible results");
    } else {
        if (lower >= computedRows || upper > computedRows)
            Rcpp::stop("bounds cannot exceed the maximum number of possible results");
    }
    
    if (IsCount) {
        if (IsGmp) {
            std::size_t sizeNum, size = sizeof(int);
            std::size_t numb = 8 * sizeof(int);
            sizeNum = sizeof(int) * (2 + (mpz_sizeinbase(computedRowMpz, 2) + numb - 1) / numb);
            size += sizeNum;
            
            SEXP ansPos = PROTECT(Rf_allocVector(RAWSXP, size));
            char* rPos = (char*)(RAW(ansPos));
            ((int*)(rPos))[0] = 1; // first int is vector-size-header
            
            // current position in rPos[] (starting after vector-size-header)
            std::size_t posPos = sizeof(int);
            posPos += myRaw(&rPos[posPos], computedRowMpz, sizeNum);
            
            Rf_setAttrib(ansPos, R_ClassSymbol, Rf_mkString("bigz"));
            UNPROTECT(1);
            return(ansPos);
        } else {
            if (computedRows > std::numeric_limits<int>::max())
                return Rcpp::wrap(computedRows);
            else
                return Rcpp::wrap(static_cast<int>(computedRows));
        }
    }
    
    std::vector<int> startZ(m);
    bool permNonTriv = false;
    double dblLower = lower;
    if (!IsGmp) mpz_set_d(lowerMpz[0], dblLower);
    
    if (bLower && mpz_cmp_ui(lowerMpz[0], 0) > 0) {
        if (IsComb) {
            if (IsGmp)
                startZ = nthCombinationGmp(n, m, lowerMpz[0], IsRepetition, IsMultiset, myReps);
            else
                startZ = nthCombination(n, m, lower, IsRepetition, IsMultiset, myReps);
        } else {
            permNonTriv = true;
            if (IsGmp)
                startZ = nthPermutationGmp(n, m, lowerMpz[0], IsRepetition, IsMultiset, myReps, freqsExpanded, true);
            else
                startZ = nthPermutation(n, m, lower, IsRepetition, IsMultiset, myReps, freqsExpanded, true);
        }
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
    
    if (IsGmp) {
        mpz_t testBound;
        mpz_init(testBound);
        
        if (bLower && bUpper) {
            mpz_sub(testBound, upperMpz[0], lowerMpz[0]);
            mpz_t absTestBound;
            mpz_init(absTestBound);
            mpz_abs(absTestBound, testBound);
            
            if (mpz_cmp_ui(absTestBound, std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(testBound);
            mpz_clear(absTestBound);
        } else if (bUpper) {
            permNonTriv = true;
            
            if (mpz_cmp_d(upperMpz[0], std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
                
            userNumRows = mpz_get_d(upperMpz[0]);
        } else if (bLower) {
            mpz_sub(testBound, computedRowMpz, lowerMpz[0]);
            mpz_abs(testBound, testBound);
            
            if (mpz_cmp_d(testBound, std::numeric_limits<int>::max()) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(testBound);
        }
        mpz_clear(testBound);
    } else {
        if (bLower && bUpper)
            userNumRows = upper - lower;
        else if (bUpper)
            userNumRows = upper;
        else if (bLower)
            userNumRows = computedRows - lower;
    }
    
    if (userNumRows == 0) {
        if (bLower && bUpper) {
            // Since lower is decremented and upper isn't, this implies that
            // upper - lower = 0 means that lower is one larger than upper as put in by the user
            
            Rcpp::stop("The number of rows must be positive. Either the lowerBound "
                           "exceeds the maximum number of possible results or the "
                           "lowerBound is greater than the upperBound.");
        } else {
            if (computedRows > std::numeric_limits<int>::max() && !IsConstrained)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = computedRows;
            
            if (!IsConstrained)
                nRows = static_cast<int>(computedRows);
        }
    } else if (userNumRows < 0) {
        Rcpp::stop("The number of rows must be positive. Either the lowerBound "
                  "exceeds the maximum number of possible results or the "
                  "lowerBound is greater than the upperBound.");
    } else if (userNumRows > std::numeric_limits<int>::max()) {
        Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
    } else {
        nRows = static_cast<int>(userNumRows);
    }
    
    const std::size_t uM = m;
    int nThreads = 1;
    
    // Determined empirically. Setting up threads can be expensive,
    // so we set the cutoff below to ensure threads aren't spawned
    // unnecessarily. We also protect users with fewer than 2 threads
    if ((nRows < 20000) || (maxThreads < 2)) {
        Parallel = false;
    } else if (!Rf_isNull(RNumThreads)) {
        int userThreads = 1;
        
        if (!Rf_isNull(RNumThreads))
            CleanConvert::convertPrimitive(RNumThreads, userThreads, "nThreads");
        
        if (userThreads > maxThreads)
            userThreads = maxThreads;
        
        // Ensure that each thread has at least 10000
        if ((nRows / userThreads) < 10000)
            userThreads = nRows / 10000;

        if (userThreads > 1 && !IsCharacter) {
            Parallel = true;
            nThreads = userThreads;
        } else {
            Parallel = false;
        }
    } else if (Parallel) {
        nThreads = (maxThreads > 2) ? (maxThreads - 1) : 2;
        
        // Ensure that each thread has at least 10000
        if ((nRows / nThreads) < 10000)
            nThreads = nRows / 10000;
    }
    
    if (IsConstrained) {
        std::vector<double> targetVals;      // numOnly = true, checkWhole = false, negPoss = true
        CleanConvert::convertVector(Rtarget, targetVals, "limitConstraints", true, false, true);
        
        if (targetVals.size() > 2)
            Rcpp::stop("there cannot be more than 2 limitConstraints");
        else if (targetVals.size() == 2 && targetVals[0] == targetVals[1])
            Rcpp::stop("The limitConstraints must be different");
        
        std::string mainFun = Rcpp::as<std::string>(f1);
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
            
            if (itComp == compForms.end())
                Rcpp::stop("comparison operators must be one of the following: '>', '>=', '<', '<=', or '=='");
                
            int myIndex = std::distance(compForms.cbegin(), itComp);
            
            // The first 5 are "standard" whereas the 6th and 7th
            // are written with the equality first. Converting
            // them here makes it easier to deal with later.
            if (myIndex > 4)
                myIndex -= 3;
            
            compFunVec[i] = compForms[myIndex];
        }
        
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
                
                bool sortNeeded = false;
                
                if (compFunVec[0].substr(0, 1) == ">" && std::min(targetVals[0], targetVals[1]) == targetVals[0]) {
                    compFunVec[0] = compFunVec[0] + "," + compFunVec[1];
                    sortNeeded = true;
                } else if (compFunVec[0].substr(0, 1) == "<" && std::max(targetVals[0], targetVals[1]) == targetVals[0]) {
                    compFunVec[0] = compFunVec[1] + "," + compFunVec[0];
                    sortNeeded = true;
                }
                
                if (sortNeeded) {
                    if (std::max(targetVals[0], targetVals[1]) == targetVals[1]) {
                        double temp = targetVals[0];
                        targetVals[0] = targetVals[1];
                        targetVals[1] = temp;
                    }
                    
                    compFunVec.pop_back();
                }
            }
        } else {
            if (targetVals.size() == 2)
                targetVals.pop_back();
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
        
        Rcpp::XPtr<funcPtr<double>> xpFunDbl = putFunPtrInXPtr<double>(mainFun);
        funcPtr<double> myFunDbl = *xpFunDbl;
        
        IsInteger = (IsInteger) && checkIsInteger(mainFun, uM, n, vNum, targetVals, myFunDbl, true);
        // Must be defined inside IsInteger check as targetVals could be outside integer data type range
        std::vector<int> targetIntVals;
        double tolerance = 0;
        
        if (IsInteger) {
            targetIntVals.assign(targetVals.cbegin(), targetVals.cend());
        } else {
            CleanConvert::convertPrimitive(Rtolerance, tolerance, "tolerance", true, false, false, true);
        }
        
        if (SpecialCase) {
            if (IsInteger) {
                return SpecCaseRet<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, nRows, keepRes, startZ, lower,
                                                        mainFun, IsMultiset, computedRows, compFunVec, targetIntVals, IsComb,
                                                        myReps, freqsExpanded, bLower, permNonTriv, userNumRows, tolerance);
            } else {
                return SpecCaseRet<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes, startZ, lower,
                                                        mainFun, IsMultiset, computedRows, compFunVec, targetVals, IsComb,
                                                        myReps, freqsExpanded, bLower, permNonTriv, userNumRows, tolerance);
            }
        }
        
        bool bUserRows = bLower || bUpper;
        
        if (mainFun == "sum" && compFunVec[0] == "==" && !IsMultiset && n > 1 && m > 1) {
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
                            return Partitions::GeneralPartitions<Rcpp::IntegerMatrix>(n, m, v64, target64, IsRepetition,
                                                                                      userNumRows, IsComb, keepRes, bUserRows);
                        } else {
                            return Partitions::GeneralPartitions<Rcpp::NumericMatrix>(n, m, v64, target64, IsRepetition,
                                                                                      userNumRows, IsComb, keepRes, bUserRows);
                        }
                    }
                }
            }
        }
        
        if (IsInteger) {
            return CombinatoricsConstraints<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, mainFun, compFunVec, targetIntVals,
                                                                 userNumRows, IsComb, keepRes, myReps, IsMultiset, tolerance, bUserRows);
        }
        
        return CombinatoricsConstraints<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, mainFun, compFunVec, targetVals,
                                                             userNumRows, IsComb, keepRes, myReps, IsMultiset, tolerance, bUserRows);
    } else {
        bool applyFun = !Rf_isNull(stdFun) && !IsFactor;

        if (applyFun) {
            if (!Rf_isFunction(stdFun))
                Rcpp::stop("FUN must be a function!");
            
            SEXP ans = PROTECT(Rf_allocVector(VECSXP, nRows));
            SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));
            
            if (IsCharacter) {
                ApplyFunction(n, m, rcppChar, IsRepetition, nRows, IsComb, myReps, 
                              ans, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            } else if (IsLogical || IsInteger) {
                Rcpp::IntegerVector rcppVInt(vInt.cbegin(), vInt.cend());
                ApplyFunction(n, m, rcppVInt, IsRepetition, nRows, IsComb, myReps, 
                              ans, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            } else {
                Rcpp::NumericVector rcppVNum(vNum.cbegin(), vNum.cend());
                ApplyFunction(n, m, rcppVNum, IsRepetition, nRows, IsComb, myReps, 
                              ans, freqsExpanded, startZ, IsMultiset, sexpFun, myEnv, 0);
            }
            
            UNPROTECT(2);
            return ans;
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
        
        std::string mainFun;
        funcPtr<double> myFunDbl;
        funcPtr<int> myFunInt;
        int nCol = m;
        
        if (keepRes) {
            mainFun = Rcpp::as<std::string>(f1);
            if (mainFun != "prod" && mainFun != "sum" && mainFun != "mean"
                    && mainFun != "max" && mainFun != "min") {
                Rcpp::stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
            }
            
            Rcpp::XPtr<funcPtr<double>> xpFunDbl = putFunPtrInXPtr<double>(mainFun);
            Rcpp::XPtr<funcPtr<int>> xpFunInt = putFunPtrInXPtr<int>(mainFun);
            myFunDbl = *xpFunDbl;
            myFunInt = *xpFunInt;
            
            IsInteger = (IsInteger) && checkIsInteger(mainFun, uM, n, vNum, vNum, myFunDbl);
            ++nCol;
        }
        
        if (Parallel) {
            permNonTriv = true;
            RcppThread::ThreadPool pool(nThreads);
            int step = 0, stepSize = nRows / nThreads;
            int nextStep = stepSize;
            
            if (IsLogical) {
                Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(nRows, nCol);
                RcppParallel::RMatrix<int> parBool(matBool);

                for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                    pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<int>, int>), 
                              n, m, vInt, IsRepetition, nextStep, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunInt, keepRes, std::ref(parBool), step);
                    
                    getStartZ(n, m, lower, stepSize, lowerMpz[0], IsRepetition,
                              IsComb, IsMultiset, IsGmp, myReps, freqsExpanded, startZ);
                }

                pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<int>, int>), 
                          n, m, vInt, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                          startZ, permNonTriv, IsMultiset, myFunInt, keepRes, std::ref(parBool), step);

                pool.join();
                return matBool;
                
            } else if (IsFactor || IsInteger) {
                Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(nRows, nCol);
                RcppParallel::RMatrix<int> parInt(matInt);
                
                for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                    pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<int>, int>), 
                              n, m, vInt, IsRepetition, nextStep, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunInt, keepRes, std::ref(parInt), step);

                    getStartZ(n, m, lower, stepSize, lowerMpz[0], IsRepetition, 
                               IsComb, IsMultiset, IsGmp, myReps, freqsExpanded, startZ);
                }
                
                pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<int>, int>),
                          n, m, vInt, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                          startZ, permNonTriv, IsMultiset, myFunInt, keepRes, std::ref(parInt), step);
                
                pool.join();
                
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
                RcppParallel::RMatrix<double> parNum(matNum);
                
                for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                    pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<double>, double>), 
                              n, m, vNum, IsRepetition, nextStep, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunDbl, keepRes, std::ref(parNum), step);
                    
                    getStartZ(n, m, lower, stepSize, lowerMpz[0], IsRepetition, 
                              IsComb, IsMultiset, IsGmp, myReps, freqsExpanded, startZ);
                }

                pool.push(std::cref(GeneralReturn<RcppParallel::RMatrix<double>, double>),
                          n, m, vNum, IsRepetition, nRows, IsComb, myReps, freqsExpanded, startZ,
                          permNonTriv, IsMultiset, myFunDbl, keepRes, std::ref(parNum), step);
                
                pool.join();
                return matNum;
            }
        } else {
            if (IsCharacter) {
                Rcpp::CharacterMatrix matChar = Rcpp::no_init_matrix(nRows, nCol);
                CharacterReturn(n, m, rcppChar, IsRepetition, nRows, IsComb, myReps,
                                freqsExpanded, startZ, permNonTriv, IsMultiset, keepRes, matChar, 0);
                return matChar;
                
            } else if (IsLogical) {
                Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(nRows, nCol);
                GeneralReturn(n, m, vInt, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunInt, keepRes, matBool, 0);
                return matBool;
                
            } else if (IsFactor || IsInteger) {
                Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(nRows, nCol);
                GeneralReturn(n, m, vInt, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunInt, keepRes, matInt, 0);
                
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
                GeneralReturn(n, m, vNum, IsRepetition, nRows, IsComb, myReps, freqsExpanded, 
                              startZ, permNonTriv, IsMultiset, myFunDbl, keepRes, matNum, 0);
                return matNum;
            }
        }
    }
}
