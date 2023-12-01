#include "cpp11/strings.hpp"

#include "Combinations/NthCombination.h"
#include "Permutations/NthPermutation.h"
#include "Permutations/PermuteCount.h"
#include "CppConvert.h"
#include <algorithm> // std:: sort, std::minmax
#include <numeric>   // std::iota
#include <cmath>

static gmp_randstate_t seed_state;
static int seed_init = 0;

void SetType(VecType &myType, SEXP Rv) {

    if (Rf_isMatrix(Rv)) {
        cpp11::stop("Matrices are not supported for v");
    }

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
            // Vectors of class bigZ and mpfr cause a lot of headaches, and for
            // this we simply exclude all raw vectors that have any attributes.
            // If you think there is a clean solution for including these
            // cases, please contact me @ jwood000@gmail.com. N.B., see commit
            // 655 which includes a function for returning a matrix of class
            // bigz. I observed terrible performance compared to simply
            // converting to a character vector.
            if (ATTRIB(Rv) == R_NilValue) {
                myType = VecType::Raw;
                break;
            }
        } default: {
            cpp11::stop("Only atomic types are supported for v");
        }
    }
}

void SetFactorClass(SEXP res, SEXP Rv) {
    Rf_setAttrib(res, R_ClassSymbol, Rf_getAttrib(Rv, R_ClassSymbol));
    Rf_setAttrib(res, R_LevelsSymbol, Rf_getAttrib(Rv, R_LevelsSymbol));
}

bool IsDecimal(SEXP Rv) {

    if (TYPEOF(Rv) == REALSXP && Rf_length(Rv) == 1) {
        double res = Rf_asReal(Rv);

        if (static_cast<std::int64_t>(res) == res) {
            return false;
        } else {
            return true;
        }
    } else {
        return false;
    }
}

int GetLength(SEXP Rv, VecType myType) {

    if (myType >= VecType::Logical) {
        return(Rf_length(Rv));
    } else if (IsDecimal(Rv)) {
        return(1);
    } else if (Rf_length(Rv) == 1) {
        int seqEnd = 0;

        // numOnly = true, checkWhole = true, negPoss = true
        CppConvert::convertPrimitive(Rv, seqEnd, myType,
                                       "v, if v is not a character"
                                       " and of length 1,",
                                       true, true, true);

        const int first = (seqEnd < 0) ? -1 : ((seqEnd == 0) ? 0 : 1);
        std::pair<int, int> mnmx = std::minmax(first, seqEnd);
        const int n = mnmx.second - mnmx.first + 1;

        if (n >= maxVecSize) {
            cpp11::stop("Not enough memory! The vector you have"
                            " requested is larger than %s",
                            std::to_string(maxVecSize).c_str());
        }

        return(n);
    } else {
        return(Rf_length(Rv));
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
        CppConvert::convertVector(RFreqs, Reps, VecType::Integer, "freqs");
        const bool allOne = std::all_of(Reps.cbegin(), Reps.cend(),
                                        [](int v_i) {return v_i == 1;});
        if (allOne) {
            IsMult = false;
            freqs.assign(Reps.size(), 1);
            Reps.clear();
        } else {
            IsMult = true;

            for (std::size_t i = 0; i < Reps.size(); ++i) {
                for (int j = 0; j < Reps[i]; ++j) {
                    freqs.push_back(i);
                }
            }
        }
    }

    if (Rf_isNull(Rm)) {
        m = freqs.empty() ? n : freqs.size();
    } else if (Rf_length(Rm) > 1) {
        cpp11::stop("length of m must be 1");
    } else {
        CppConvert::convertPrimitive(Rm, m, VecType::Integer, "m");
    }
}

// This function is mainly for handling missing data. When a
// constraint is applied, we must throw out these values. Also
// note that the constraint algos are expecting sorted data
void SetFinalValues(VecType &myType, std::vector<int> &Reps,
                    std::vector<int> &freqs, std::vector<int> &vInt,
                    std::vector<double> &vNum, int &n, int &m,
                    bool IsMult, bool IsRep, bool IsConstrained) {

    if (IsConstrained && vNum.size()) {
        bool any_na = false;

        for (int i = (vNum.size() - 1); i >= 0; --i) {
            if (CppConvert::CheckNA(vNum[i], myType)) {
                any_na = true;
                vNum.erase(vNum.begin() + i);

                if (IsMult) {
                    Reps.erase(Reps.begin() + i);
                }
            }
        }

        if (IsRep) {
            vNum.erase(std::unique(vNum.begin(), vNum.end()), vNum.end());
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

        if (any_na) {
            freqs.clear();

            for (int i = 0; i < static_cast<int>(Reps.size()); ++i) {
                for (int j = 0; j < Reps[i]; ++j) {
                    freqs.push_back(i);
                }
            }
        }
    }

    if (myType == VecType::Integer && vInt.size() == 0) {
        vInt.clear();
        vInt.reserve(n);

        for (auto v_i: vNum) {
            vInt.push_back(static_cast<int>(v_i));
        }
    }

    if (IsMult) {
        // The 'freqs' in the error message below refers to the user
        // supplied parameter named 'freqs'. When we use freqs as a
        // C++ variable, we are referring to the expansion of this
        // variable. myReps and Reps in C++ are converted directly
        // from the user supplied 'freqs' parameter.
        if (n != static_cast<int>(Reps.size())) {
            cpp11::stop("the length of freqs must equal the length of v");
        }

        if (m > static_cast<int>(freqs.size())) {
            m = freqs.size();
        }
    } else if (!IsRep && m > n) {
        cpp11::stop("m must be less than or equal to the length of v");
    }
}

void SetBasic(SEXP Rv, std::vector<double> &vNum,
              std::vector<int> &vInt, int &n, VecType &myType) {

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
        CppConvert::convertPrimitive(
            Rv, seqEnd, myType, "v, if v is not a character"
            " and of length 1,", true, true, true
        );

        const int first = (seqEnd < 0) ? -1 : ((seqEnd == 0) ? 0 : 1);
        std::pair<int, int> mnmx = std::minmax(first, seqEnd);
        n = mnmx.second - mnmx.first + 1;

        if (n >= maxVecSize) {
            cpp11::stop("Not enough memory! The vector you have"
                        " requested is larger than %s",
                        std::to_string(maxVecSize).c_str());
        }

        vNum.resize(n);
        std::iota(vNum.begin(), vNum.end(), static_cast<double>(mnmx.first));
    } else {
        vNum = CppConvert::GetNumVec<double>(Rv);
        n = vNum.size();
    }
}

void SetValues(VecType &myType, std::vector<int> &Reps,
               std::vector<int> &freqs, std::vector<int> &vInt,
               std::vector<double> &vNum, SEXP Rv, SEXP RFreqs,
               SEXP Rm, int &n, int &m, bool &IsMult,
               bool &IsRep, bool IsConstrained) {

    SetBasic(Rv, vNum, vInt, n, myType);
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

        if (!Rf_isNull(RNumThreads)) {
            CppConvert::convertPrimitive(RNumThreads, userThreads,
                                         VecType::Integer, "nThreads");
        }

        if (userThreads > maxThreads) {
            userThreads = maxThreads;
        }

        // Ensure that each thread has at least halfLimit
        if ((nRows / userThreads) < halfLimit) {
            userThreads = nRows / halfLimit;
        }

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
        if ((nRows / nThreads) < halfLimit) {
            nThreads = nRows / halfLimit;
        }
    }
}

void SetNumResults(bool IsGmp, bool bLower, bool bUpper, bool bSetNum,
                   const mpz_class &upperMpz, const mpz_class &lowerMpz,
                   double lower, double upper, double computedRows,
                   const mpz_class &computedRowsMpz, int &nRows,
                   double &userNumRows) {

    if (IsGmp) {
        mpz_class testBound;

        if (bLower && bUpper) {
            testBound = upperMpz - lowerMpz;
            mpz_class absTestBound(abs(testBound));

            if (cmp(absTestBound, std::numeric_limits<int>::max()) > 0) {
                cpp11::stop("The number of rows cannot exceed 2^31 - 1.");
            }

            userNumRows = testBound.get_d();
        } else if (bUpper) {
            if (cmp(upperMpz, std::numeric_limits<int>::max()) > 0) {
                cpp11::stop("The number of rows cannot exceed 2^31 - 1.");
            }

            userNumRows = upperMpz.get_d();
        } else if (bLower) {
            testBound = abs(computedRowsMpz - lowerMpz);

            if (cmp(testBound, std::numeric_limits<int>::max()) > 0) {
                cpp11::stop("The number of rows cannot exceed 2^31 - 1.");
            }

            userNumRows = testBound.get_d();
        }
    } else {
        if (bLower && bUpper) {
            userNumRows = upper - lower;
        } else if (bUpper) {
            userNumRows = upper;
        } else if (bLower) {
            userNumRows = computedRows - lower;
        }
    }

    if (userNumRows == 0) {
        if (bLower && bUpper) {
            // Since lower is decremented and upper isn't, this implies that upper - lower = 0
            // which means that lower is one larger than upper as put in by the user

            cpp11::stop("The number of rows must be positive. Either the"
                     "lowerBound exceeds the maximum number of possible"
                     " results or the lowerBound is greater "
                     "than the upperBound.");
        } else {
            // See comment in ConstraintsMain.cpp. Basically, we don't want to
            // throw an error when we don't really know how many constrained
            // results we have as computedRows is a strict upper bound and not
            // the least upper bound.
            if (bSetNum && computedRows > std::numeric_limits<int>::max()) {
                cpp11::stop("The number of rows cannot exceed 2^31 - 1.");
            }

            userNumRows = computedRows;

            if (bSetNum) {
                nRows = static_cast<int>(computedRows);
            }
        }
    } else if (userNumRows < 0) {
        cpp11::stop("The number of rows must be positive. Either the lowerBound"
                 " exceeds the maximum number of possible results or the"
                 " lowerBound is greater than the upperBound.");
    } else if (userNumRows > std::numeric_limits<int>::max()) {
        cpp11::stop("The number of rows cannot exceed 2^31 - 1.");
    } else {
        nRows = static_cast<int>(userNumRows);
    }
}

void SetBounds(SEXP Rlow, SEXP Rhigh, bool IsGmp, bool &bLower,
               bool &bUpper, double &lower, double &upper,
               mpz_class &lowerMpz, mpz_class &upperMpz,
               const mpz_class &computedRowsMpz, double computedRows) {

    if (!Rf_isNull(Rlow)) {
        if (IsGmp) {
            CppConvert::convertMpzClass(Rlow, lowerMpz, "lower");
            bLower = cmp(lowerMpz, 1) > 0;
            lower = bLower ? 1 : 0;

            if (cmp(lowerMpz, computedRowsMpz) > 0) {
                cpp11::stop("bounds cannot exceed the maximum "
                             "number of possible results");
            }

            --lowerMpz;
        } else {                                    // numOnly = false
            CppConvert::convertPrimitive(Rlow, lower, VecType::Numeric,
                                         "lower", false);
            bLower = lower > 1;

            if (lower > computedRows) {
                cpp11::stop("bounds cannot exceed the maximum "
                             "number of possible results");
            }

            --lower;
        }
    }

    if (!Rf_isNull(Rhigh)) {
        bUpper = true;

        if (IsGmp) {
            CppConvert::convertMpzClass(Rhigh, upperMpz, "upper");

            if (cmp(upperMpz, computedRowsMpz) > 0) {
                cpp11::stop("bounds cannot exceed the maximum "
                             "number of possible results");
            }
        } else {
            CppConvert::convertPrimitive(Rhigh, upper,   // numOnly = false
                                         VecType::Numeric, "upper", false);

            if (upper > computedRows) {
                cpp11::stop("bounds cannot exceed the maximum "
                             "number of possible results");
            }
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

                if (phaseOneDbl * (static_cast<double>(m) - first) >
                        sizeTestVec.max_size()) {
                    generalRet = true;
                }
            }
        }

        phaseOne = (generalRet) ? 0 : static_cast<int>(phaseOneDbl);
    }
}

void SetStartZ(const std::vector<int> &myReps,
               const std::vector<int> &freqs, std::vector<int> &z,
               bool IsComb, int n, int m, double lower,
               const mpz_class &lowerMpz, bool IsRep,
               bool IsMult, bool IsGmp) {

    if (lower > 0) {
        if (IsComb) {
            const nthCombPtr nthCombFun = GetNthCombFunc(IsMult,
                                                         IsRep, IsGmp);
            z = nthCombFun(n, m, lower, lowerMpz, myReps);
        } else {
            const nthPermPtr nthPermFun = GetNthPermFunc(IsMult,
                                                         IsRep, IsGmp);
            z = nthPermFun(n, m, lower, lowerMpz, myReps);
            TopOffPerm(z, myReps, n, m, IsRep, IsMult);
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

void SetRandomSample(SEXP RindexVec, SEXP RNumSamp, int &sampSize,
                     bool IsGmp, double computedRows,
                     std::vector<double> &mySample,
                     SEXP baseSample, SEXP rho) {

    // We must treat gmp case special. We first have to get the size of our sample
    // vector, as we have to declare a mpz_class array with known size. You will note
    // that in the base case below, we simply populate mySample, otherwise we just
    // get the size. This size var will be used in the next block (If (IsGmp)...)
    if (Rf_isNull(RindexVec)) {
        if (Rf_isNull(RNumSamp)) {
            cpp11::stop("n and sampleVec cannot both be NULL");
        }

        if (Rf_length(RNumSamp) > 1) {
            cpp11::stop("length of n must be 1. For specific "
                        "combinations, use sampleVec.");
        }

        CppConvert::convertPrimitive(RNumSamp, sampSize,
                                       VecType::Integer, "n");

        if (!IsGmp) {
            if (sampSize > computedRows) {
                cpp11::stop("n exceeds the maximum number of possible results");
            }

            cpp11::sexp sample = Rf_lang3(baseSample,
                                          Rf_ScalarReal(computedRows),
                                          RNumSamp);
            cpp11::sexp val = Rf_eval(sample, rho);
            mySample.resize(sampSize);

            if (computedRows < std::numeric_limits<int>::max()) {
                int* intSamp = INTEGER(val);

                for (int j = 0; j < sampSize; ++j) {
                    mySample[j] = intSamp[j];
                }
            } else {
                double* dblSamp = REAL(val);

                for (int j = 0; j < sampSize; ++j) {
                    mySample[j] = dblSamp[j];
                }
            }

        }
    } else {
        if (IsGmp) {
            switch (TYPEOF(RindexVec)) {
                case RAWSXP: {
                    const char* raw = (char*) RAW(RindexVec);
                    sampSize = ((int*) raw)[0];
                    break;
                } default: {
                    sampSize = LENGTH(RindexVec);
                }
            }
        } else {
            CppConvert::convertVector(RindexVec, mySample,
                                        VecType::Numeric,
                                        "sampleVec", false);
            sampSize = mySample.size();
            const double myMax = *std::max_element(mySample.cbegin(),
                                                   mySample.cend());

            if (myMax > computedRows) {
                cpp11::stop("One or more of the requested values in sampleVec "
                         "exceeds the maximum number of possible results");
            }
        }

        if (sampSize > std::numeric_limits<int>::max()) {
            cpp11::stop("The number of rows cannot exceed 2^31 - 1");
        }
    }

    // Get zero base index
    for (auto &s: mySample) {
        --s;
    }
}

void SetRandomSampleMpz(SEXP RindexVec, SEXP RmySeed, int sampSize,
                        bool IsGmp, const mpz_class &computedRowsMpz,
                        std::vector<mpz_class> &myVec) {

    if (IsGmp) {
        if (!Rf_isNull(RindexVec)) {
            CppConvert::convertMPZVector(RindexVec, myVec,
                                         sampSize, "sampleVec");

            // get zero base
            for (int i = 0; i < sampSize; ++i) {
                --myVec[i];
            }
        } else {
            // The following code is very similar to the source
            // code of gmp::urand.bigz. The main difference is
            // the use of mpz_urandomm instead of mpz_urandomb
            if (seed_init == 0) {
                gmp_randinit_default(seed_state);
            }

            seed_init = 1;

            if (!Rf_isNull(RmySeed)) {
                mpz_class mpzSeed;
                CppConvert::convertMpzClass(RmySeed, mpzSeed, "seed");
                gmp_randseed(seed_state, mpzSeed.get_mpz_t());
            }

            // random number is between 0 and gmpRows[0] - 1
            // so we need to add 1 to each element
            for (int i = 0; i < sampSize; ++i) {
                mpz_urandomm(myVec[i].get_mpz_t(), seed_state,
                             computedRowsMpz.get_mpz_t());
            }
        }

        mpz_class maxGmp(myVec.front());

        for (int i = 1; i < sampSize; ++i) {
            if (cmp(myVec[i], maxGmp) > 0) {
                maxGmp = myVec[i];
            }
        }

        if (cmp(maxGmp, computedRowsMpz) >= 0) {
            cpp11::stop("One or more of the requested values in sampleVec "
                     "exceeds the maximum number of possible results");
        }
    }
}

void SetSampleNames(SEXP object, bool IsGmp, int sampSize,
                    const std::vector<double> &mySample,
                    const std::vector<mpz_class> &myBigSamp,
                    bool IsNamed, SEXP colNames, int xtraDims) {

    if (IsNamed) {
        cpp11::writable::strings myNames(sampSize);

        if (IsGmp) {
            mpz_class temp;

            for (int i = 0; i < sampSize; ++i) {
                temp = myBigSamp[i] + 1;
                myNames[i] = temp.get_str();
            }
        } else {
            for (int i = 0; i < sampSize; ++i) {
                const std::string name = std::to_string(
                    static_cast<int64_t>(mySample[i] + 1)
                );

                myNames[i] = name.c_str();
            }
        }

        if (Rf_isMatrix(object) || Rf_isArray(object)) {
            cpp11::sexp dimNames = Rf_allocVector(VECSXP, 1 + xtraDims);
            SET_VECTOR_ELT(dimNames, 0, myNames);
            if (xtraDims) SET_VECTOR_ELT(dimNames, xtraDims, colNames);
            Rf_setAttrib(object, R_DimNamesSymbol, dimNames);
        } else if (Rf_isList(object) || Rf_isVector(object)) {
            Rf_setAttrib(object, R_NamesSymbol, myNames);
        }
    }
}

SEXP GetInt64Vec(const std::vector<std::int64_t> &v) {
    const int size = v.size();
    cpp11::sexp res = Rf_allocVector(REALSXP, size);
    double* ptrRes = REAL(res);

    for (int i = 0; i < size; ++i) {
        ptrRes[i] = v[i];
    }

    return res;
}
