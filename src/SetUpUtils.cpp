#include "Combinations/NthCombination.h"
#include "Permutations/NthPermutation.h"
#include "Permutations/PermuteCount.h"
#include "ImportExportMPZ.h"
#include "CleanConvert.h"
#include <algorithm>
#include <numeric>

void SetType(VecType &myType, SEXP Rv) {
    
    switch(TYPEOF(Rv)) {
        case LGLSXP: {
            myType = VecType::Logical;
            break;
        } case INTSXP: {
            myType = VecType::Integer;
            break;
        } case REALSXP: {
            myType = VecType::Numeric;
            break;
        } case STRSXP: {
            myType = VecType::Character;
            break;
        } case CPLXSXP: {
            myType = VecType::Complex;
            break;
        } case RAWSXP: {
            // Vectors of class bigZ and mpfr cause a lot of headaches, and for this
            // we simply exclude all raw vectors that have any attributes. If you
            // think there is a clean solution for including these cases, please
            // contact me @ jwood000@gmail.com. N.B., see commit 655 which includes
            // a function for returning a matrix of class bigz. I observed terrible
            // performance compared to simply converting to a character vector.
            if (ATTRIB(Rv) == R_NilValue) {
                myType = VecType::Raw;
                break;
            }
        } default: {
            Rf_error("Only atomic types are supported for v");
        }
    }
}

void SetFactorClass(SEXP res, SEXP Rv) {
    SEXP myLevels = Rf_getAttrib(Rv, R_LevelsSymbol);
    SEXP myClass = Rf_getAttrib(Rv, R_ClassSymbol);
    Rf_setAttrib(res, R_ClassSymbol, myClass);
    Rf_setAttrib(res, R_LevelsSymbol, myLevels);
}

bool IsDecimal(SEXP Rv) {
    
    if (TYPEOF(Rv) == REALSXP && Rf_length(Rv) == 1) {
        double res = REAL(Rv)[0];
        
        if (static_cast<std::int64_t>(res) == res) {
            return false;
        } else {
            return true;
        }
    } else {
        return false;
    }
}

void SetFreqsAndM(std::vector<int> &Reps,
                  std::vector<int> &freqs, SEXP RFreqs, SEXP Rm,
                  int &n, int &m, bool &IsMult, bool &IsRep) {
    
    if (Rf_isNull(RFreqs)) {
        IsMult = false;
    } else {
        // If user passes repetition = TRUE as well as non-trivial freqs, we give
        // preference to freqs as user may assume that since multisets includes
        // replication of certain elements, then repetition must be set to TRUE.
        IsRep = false;
        CleanConvert::convertVector(RFreqs, Reps, VecType::Integer, "freqs");
        const bool allOne = std::all_of(Reps.cbegin(), Reps.cend(), 
                                        [](int v_i) {return v_i == 1;});
        if (allOne) {
            IsMult = false;
            freqs.assign(Reps.size(), 1);
            Reps.clear();
        } else {
            IsMult = true;
            
            for (std::size_t i = 0; i < Reps.size(); ++i)
                for (int j = 0; j < Reps[i]; ++j)
                    freqs.push_back(i);
        }
    }
    
    if (Rf_isNull(Rm)) {
        m = (freqs.empty()) ? n : freqs.size();
    } else {
        if (Rf_length(Rm) > 1)
            Rf_error("length of m must be 1");
        
        CleanConvert::convertPrimitive(Rm, m, VecType::Integer, "m");
    }
}

// This function is mainly for handling missing data. When a
// constraint is applied, we must throw out these values. Also
// note that the constraint algos are expecting sorted data
void SetFinalValues(VecType &myType, std::vector<int> &Reps,
                    std::vector<int> &freqs, std::vector<int> &vInt,
                    std::vector<double> &vNum, int &n, int &m,
                    bool IsMult, bool IsRep, bool IsConstrained) {
    
    if (IsConstrained) {
        for (int i = (vNum.size() - 1); i >= 0; --i) {
            if (CleanConvert::CheckNA(vNum[i], myType)) {
                vNum.erase(vNum.begin() + i);
                if (IsMult) Reps.erase(Reps.begin() + i);
            }
        }
        
        n = vNum.size();
        
        if (IsMult) {
            for (int i = 0; i < (n - 1); ++i) {
                for (int j = (i + 1); j < n; ++j) {
                    if (vNum[i] > vNum[j]) {
                        std::swap(vNum[i], vNum[j]);
                        std::swap(Reps[i], Reps[j]);
                    }
                }
            }
        } else {
            std::sort(vNum.begin(), vNum.end());
        }
        
        freqs.clear();
        
        for (int i = 0; i < static_cast<int>(Reps.size()); ++i)
            for (int j = 0; j < Reps[i]; ++j)
                freqs.push_back(i);
    } else {
        VecType OldType = myType;
        
        for (int i = 0; i < n; ++i) {
            if (CleanConvert::CheckNA(vNum[i], OldType)) {
                myType = VecType::Numeric;
                
                if (OldType == VecType::Integer) {
                    vNum[i] = NA_REAL;
                }
            }
        }
    }
    
    if (myType == VecType::Integer) {
        vInt.assign(vNum.cbegin(), vNum.cend());
    }
    
    if (IsMult) {
        // The 'freqs' in the error message below refers to the user
        // supplied parameter named 'freqs'. When we use freqs as a
        // C++ variable, we are referring to the expansion of this
        // variable. myReps and Reps in C++ are converted directly
        // from the user supplied 'freqs' parameter.
        if (n != Reps.size())
            Rf_error("the length of freqs must equal the length of v");
        
        if (m > static_cast<int>(freqs.size())) {
            m = freqs.size();
        }
    } else if (!IsRep && m > n) {
        Rf_error("m must be less than or equal to the length of v");
    }
}

void SetValues(VecType &myType, std::vector<int> &Reps,
               std::vector<int> &freqs, std::vector<int> &vInt,
               std::vector<double> &vNum, SEXP Rv, SEXP RFreqs,
               SEXP Rm, int &n, int &m, bool &IsMult,
               bool &IsRep, bool IsConstrained) {
    
    if (myType > VecType::Logical) {
        n = Rf_length(Rv);
    } else if (IsDecimal(Rv)) {
        vNum.resize(1);
        vNum[0] = REAL(Rv)[0];
        n = 1;
    } else if (myType == VecType::Logical) {
        int* intVec = INTEGER(Rv);
        n = Rf_length(Rv);
        vInt.assign(intVec, intVec + n);
    } else if (Rf_length(Rv) == 1) {
        int seqEnd = 0;
        myType = VecType::Integer;
        
        // numOnly = true, checkWhole = true, negPoss = true
        CleanConvert::convertPrimitive(Rv, seqEnd, myType,
                                       "If v is not a character"
                                       " and of length 1, it",
                                       true, true, true);
        
        std::pair<int, int> mnmx = std::minmax(1, seqEnd);
        n = mnmx.second - mnmx.first + 1;
        constexpr int maxVecSize = std::numeric_limits<int>::max() / 2;
        
        if (n < maxVecSize) {
            vNum.resize(n);
        } else {
            Rf_error("Not enough memory! The vector you have"
                         " requested is larger than %s",
                         std::to_string(maxVecSize).c_str());
        }
        
        std::iota(vNum.begin(), vNum.end(), mnmx.first);
    } else {
        vNum = CleanConvert::GetNumVec<double>(Rv);
        n = vNum.size();
    }
    
    SetFreqsAndM(Reps, freqs, RFreqs, Rm, n, m, IsMult, IsRep);
    SetFinalValues(myType, Reps, freqs, vInt, vNum,
                   n, m, IsMult, IsRep, IsConstrained);
}

void SetThreads(bool &Parallel, int maxThreads, int nRows,
                VecType myType, int &nThreads, SEXP RNumThreads, int limit) {
    
    const int halfLimit = limit / 2;
    
    // Determined empirically. Setting up threads can be expensive,
    // so we set the cutoff below to ensure threads aren't spawned
    // unnecessarily. We also protect users with fewer than 2 threads
    if ((nRows < limit) || (maxThreads < 2) || myType > VecType::Logical) {
        Parallel = false;
    } else if (!Rf_isNull(RNumThreads)) {
        int userThreads = 1;
        
        if (!Rf_isNull(RNumThreads))
            CleanConvert::convertPrimitive(RNumThreads, userThreads,
                                           VecType::Integer, "nThreads");
        
        if (userThreads > maxThreads)
            userThreads = maxThreads;
        
        // Ensure that each thread has at least halfLimit
        if ((nRows / userThreads) < halfLimit)
            userThreads = nRows / halfLimit;
        
        if (userThreads > 1) {
            Parallel = true;
            nThreads = userThreads;
        } else {
            Parallel = false;
        }
    } else if (Parallel) {
        // We have already ruled out cases when the user has fewer than 2 
        // threads. So if user has exactly 2 threads, we enable them both.
        nThreads = (maxThreads > 2) ? (maxThreads - 1) : maxThreads;
        
        // Ensure that each thread has at least halfLimit
        if ((nRows / nThreads) < halfLimit)
            nThreads = nRows / halfLimit;
    }
}

void SetNumResults(bool IsGmp, bool bLower, bool bUpper, bool IsGenCnstrd,
                   mpz_t *const upperMpz, mpz_t *const lowerMpz, double lower,
                   double upper, double computedRows, mpz_t &computedRowsMpz,
                   int &nRows, double &userNumRows) {
    
    if (IsGmp) {
        mpz_t testBound;
        mpz_init(testBound);
        
        if (bLower && bUpper) {
            mpz_sub(testBound, upperMpz[0], lowerMpz[0]);
            mpz_t absTestBound;
            mpz_init(absTestBound);
            mpz_abs(absTestBound, testBound);
            
            if (mpz_cmp_ui(absTestBound, std::numeric_limits<int>::max()) > 0)
                Rf_error("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(testBound);
            mpz_clear(absTestBound);
        } else if (bUpper) {
            if (mpz_cmp_d(upperMpz[0], std::numeric_limits<int>::max()) > 0)
                Rf_error("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(upperMpz[0]);
        } else if (bLower) {
            mpz_sub(testBound, computedRowsMpz, lowerMpz[0]);
            mpz_abs(testBound, testBound);
            
            if (mpz_cmp_d(testBound, std::numeric_limits<int>::max()) > 0)
                Rf_error("The number of rows cannot exceed 2^31 - 1.");
            
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
            // Since lower is decremented and upper isn't, this implies that upper - lower = 0
            // which means that lower is one larger than upper as put in by the user
            
            Rf_error("The number of rows must be positive. Either the lowerBound "
                           "exceeds the maximum number of possible results or the "
                           "lowerBound is greater than the upperBound.");
        } else {
            if (computedRows > std::numeric_limits<int>::max() && !IsGenCnstrd)
                Rf_error("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = computedRows;
            
            if (!IsGenCnstrd)
                nRows = static_cast<int>(computedRows);
        }
    } else if (userNumRows < 0) {
        Rf_error("The number of rows must be positive. Either the lowerBound "
                       "exceeds the maximum number of possible results or the "
                       "lowerBound is greater than the upperBound.");
    } else if (userNumRows > std::numeric_limits<int>::max()) {
        Rf_error("The number of rows cannot exceed 2^31 - 1.");
    } else {
        nRows = static_cast<int>(userNumRows);
    }
}

void SetBounds(SEXP Rlow, SEXP Rhigh, bool IsGmp, bool &bLower, bool &bUpper,
               double &lower, double &upper, mpz_t *const lowerMpz, 
               mpz_t *const upperMpz, mpz_t computedRowsMpz, double computedRows) {
    
    if (!Rf_isNull(Rlow)) {
        if (IsGmp) {
            createMPZArray(Rlow, lowerMpz, 1, "lower");
            bLower = mpz_cmp_si(lowerMpz[0], 1) > 0;
            lower = (bLower) ? 1 : 0;
            
            if (mpz_cmp(lowerMpz[0], computedRowsMpz) > 0)
                Rf_error("bounds cannot exceed the maximum number of possible results");
            
            mpz_sub_ui(lowerMpz[0], lowerMpz[0], 1);
        } else {                                    // numOnly = false
            CleanConvert::convertPrimitive(Rlow, lower,
                                           VecType::Numeric, "lower", false);
            bLower = lower > 1;
            
            if (lower > computedRows)
                Rf_error("bounds cannot exceed the maximum number of possible results");
            
            --lower;
        }
    }
    
    if (!Rf_isNull(Rhigh)) {
        bUpper = true;
        
        if (IsGmp) {
            createMPZArray(Rhigh, upperMpz, 1, "upper");
            
            if (mpz_cmp(upperMpz[0], computedRowsMpz) > 0)
                Rf_error("bounds cannot exceed the maximum number of possible results");
            
        } else {                                     // numOnly = false
            CleanConvert::convertPrimitive(Rhigh, upper,
                                           VecType::Numeric, "upper", false);
            
            if (upper > computedRows)
                Rf_error("bounds cannot exceed the maximum number of possible results");
        }
    }
}

void PermuteSpecific(int &phaseOne, bool &generalRet, int n, int m,
                     int nRows, bool IsMult, bool IsCharacter,
                     bool IsComb, bool bLower, bool IsRep) {
    
    if (!IsComb) {
        double phaseOneDbl = 0.0;
        
        // If we have permutations and trivial starting index (i.e. bLower = false)
        if (!bLower) {
            if (IsRep) {
                phaseOneDbl = std::pow(static_cast<double>(n),
                                       static_cast<double>(m - 1));
            } else {
                phaseOneDbl = NumPermsNoRep(n - 1, m - 1);
            }
        }
        
        generalRet = IsMult      ||
                     bLower      || 
                     n == 1      ||
                     IsCharacter ||
                     phaseOneDbl > std::numeric_limits<int>::max();
        
        if (!generalRet) {
            // Since we create an indexing matrix, we want to make sure
            // it is worthwhile. If phaseOne takes up too much time, we
            // are better off generating permutations one at a time,
            // thus we would need generalRet = true.
            
            if ((phaseOneDbl * 2.0) > nRows) {
                generalRet = true;
            } else {
                
                // Here, we estimate the maximum size of an array of ints
                // by taking advantage of the max_size method of vectors.
                // If the indexing matrix takes up too much space, it will
                // not be worth it, thus we fall back to generating one
                // at a time (i.e. generalRet = true).
                
                std::vector<int> sizeTestVec;
                const double first = (IsRep) ? 1.0 : 0.0;
                
                if (phaseOneDbl * (static_cast<double>(m) - first) > sizeTestVec.max_size()) {
                    generalRet = true;
                }
            }
        }
        
        phaseOne = (generalRet) ? 0 : static_cast<int>(phaseOneDbl);
    }
}

void SetStartZ(const std::vector<int> &myReps,
               const std::vector<int> &freqs,
               std::vector<int> &z, bool IsComb, int n,
               int m, double lower, mpz_t lowerMpz,
               bool IsRep, bool IsMult, bool IsGmp) {
    
    if (lower > 0) {
        if (IsComb) {
            const nthCombPtr nthCombFun = GetNthCombFunc(IsMult, IsRep, IsGmp);
            z = nthCombFun(n, m, lower, lowerMpz, myReps);
        } else {
            const nthPermPtr nthPermFun = GetNthPermFunc(IsMult, IsRep, IsGmp);
            SetStartPerm(z, nthPermFun, myReps, n, m,
                         lower, lowerMpz, IsRep, IsMult);
        }
    } else {
        if (IsComb) {
            if (IsMult)
                z.assign(freqs.cbegin(), freqs.cbegin() + m);
            else if (IsRep)
                std::fill(z.begin(), z.end(), 0);
            else
                std::iota(z.begin(), z.end(), 0);
        } else {
            if (IsMult) {
                z = freqs;
            } else if (IsRep) {
                std::fill(z.begin(), z.end(), 0);
            } else {
                z.resize(n);
                std::iota(z.begin(), z.end(), 0);
            }
        }
    }
}