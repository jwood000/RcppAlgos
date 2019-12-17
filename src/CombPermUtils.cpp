#include "GmpCombPermUtils.h"

void SetClass(VecType &myType, const SEXP &Rv) {
    
    switch(TYPEOF(Rv)) {
        case LGLSXP: {
            myType = VecType::Logical;
            break;
        }
        case INTSXP: {
            myType = VecType::Integer;
            break;
        }
        case REALSXP: {
            myType = VecType::Numeric;
            break;
        }
        case STRSXP: {
            myType = VecType::Character;
            break;
        }
        case CPLXSXP: {
            myType = VecType::Complex;
            break;
        }
        case RAWSXP: {
            // Vectors of class bigZ and mpfr cause a lot of headaches, and for this
            // we simply exclude all raw vectors that have any attributes. If you
            // think there is a clean solution for including these cases, please
            // contact me @ jwood000@gmail.com. N.B., see commit 655 which includes
            // a function for returning a matrix of class bigz. I observed terrible
            // performance compared to simply converting to a character vector.
            Rcpp::RawVector tempRaw(Rcpp::clone(Rv));
            
            if (tempRaw.attributeNames().size() == 0) {
                myType = VecType::Raw;
                break;
            }
        }
        default: {
            Rcpp::stop("Only atomic types are supported for v");   
        }
    }
}

void SetFactorClass(Rcpp::IntegerMatrix &matInt, const SEXP &Rv) {
    Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
    Rcpp::CharacterVector myClass = testFactor.attr("class");
    Rcpp::CharacterVector myLevels = testFactor.attr("levels");
    matInt.attr("class") = myClass;
    matInt.attr("levels") = myLevels;
}

void SetValues(VecType &myType, std::vector<int> &vInt,
               std::vector<double> &vNum, int &n, const SEXP &Rv) {
    
    if (myType > VecType::Logical) {
        n = Rf_length(Rv);
    } else if (myType == VecType::Logical) {
        vInt = Rcpp::as<std::vector<int>>(Rv);
        n = vInt.size();
    } else {
        if (Rf_length(Rv) == 1) {
            int seqEnd, m1, m2;         // numOnly = true, checkWhole = true, negPoss = true
            CleanConvert::convertPrimitive(Rv, seqEnd, "If v is not a character and of length 1, it", true, true, true);
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            Rcpp::IntegerVector vTemp = Rcpp::seq(m1, m2);
            vNum = Rcpp::as<std::vector<double>>(vTemp);
            myType = VecType::Integer;
        } else {
            vNum = Rcpp::as<std::vector<double>>(Rv);
        }
        
        n = vNum.size();
    }
    
    if (myType == VecType::Integer) {
        bool IsInteger = true;
        
        for (int i = 0; i < n && IsInteger; ++i) {
            if (Rcpp::NumericVector::is_na(vNum[i])) {
                IsInteger = false;
                myType = VecType::Numeric;
            }
        }
        
        if (IsInteger)
            vInt.assign(vNum.cbegin(), vNum.cend());
    }
}

void SetFreqsAndM(SEXP RFreqs, bool &IsMultiset, std::vector<int> &Reps, bool &IsRepetition,
                  int &lenFreqs, std::vector<int> &freqsExpanded, const SEXP &Rm, int n, int &m) {
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        Reps.push_back(1);
    } else {
        IsRepetition = false;
        CleanConvert::convertVector(RFreqs, Reps, "freqs");
        lenFreqs = static_cast<int>(Reps.size());
        bool allOne = std::all_of(Reps.cbegin(), Reps.cend(), 
                                  [](int v_i) {return v_i == 1;});
        if (allOne) {
            IsMultiset = false;
            freqsExpanded = Reps;
        } else {
            IsMultiset = true;
            
            for (int i = 0; i < lenFreqs; ++i)
                for (int j = 0; j < Reps[i]; ++j)
                    freqsExpanded.push_back(i);
        }
    }
    
    if (Rf_isNull(Rm)) {
        if (freqsExpanded.empty()) {
            m = n;
        } else {
            m = freqsExpanded.size();
        }
    } else {
        if (Rf_length(Rm) > 1)
            Rcpp::stop("length of m must be 1");
        
        CleanConvert::convertPrimitive(Rm, m, "m");
    }
}

SEXP CopyRv(const SEXP &Rv, const std::vector<int> &vInt,
            const std::vector<double> &vNum, VecType myType) {
    
    if (myType > VecType::Numeric) {
        return Rcpp::clone(Rv);
    } else if (myType == VecType::Integer) {
        return Rcpp::wrap(vInt);
    } else {
        return Rcpp::wrap(vNum);
    }
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
            CleanConvert::convertPrimitive(RNumThreads, userThreads, "nThreads");
        
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

void SetRandomSample(SEXP RindexVec, SEXP RNumSamp, std::size_t &sampSize,
                     bool IsGmp, double computedRows, std::vector<double> &mySample,
                     Rcpp::Function baseSample) {
    
    // We must treat gmp case special. We first have to get the size of
    // our sample vector, as we have to declare a mpz_t array with
    // known size. You will note that in the base case below,
    // we simply populate mySample, otherwise we just get the size.
    // This size var will be used in the next block... If (IsGmp)
    
    if (Rf_isNull(RindexVec)) {
        if (Rf_isNull(RNumSamp))
            Rcpp::stop("n and sampleVec cannot both be NULL");
        
        if (Rf_length(RNumSamp) > 1)
            Rcpp::stop("length of n must be 1. For specific combinations, use sampleVec.");
        
        int nPass;
        CleanConvert::convertPrimitive(RNumSamp, nPass, "n");
        sampSize = static_cast<std::size_t>(nPass);
        
        if (!IsGmp) {
            if (nPass > computedRows)
                Rcpp::stop("n exceeds the maximum number of possible results");
            
            Rcpp::NumericVector tempSamp = baseSample(computedRows, nPass);
            mySample = Rcpp::as<std::vector<double>>(tempSamp);
        }
    } else {
        if (IsGmp) {
            switch (TYPEOF(RindexVec)) {
                case RAWSXP: {
                    const char* raw = (char*) RAW(RindexVec);
                    sampSize = ((int*) raw)[0];
                    break;
                }
                default: {
                    sampSize = LENGTH(RindexVec);
                }
            }
        } else {                                             // numOnly = false
            CleanConvert::convertVector(RindexVec, mySample, "sampleVec", false);
            sampSize = mySample.size();
            
            double myMax = *std::max_element(mySample.cbegin(), mySample.cend());
            
            if (myMax > computedRows) {
                Rcpp::stop("One or more of the requested values in sampleVec "
                               "exceeds the maximum number of possible results");
            }
        }
        
        if (sampSize > std::numeric_limits<int>::max())
            Rcpp::stop("The number of rows cannot exceed 2^31 - 1");
    }
    
    // Get zero base index
    for (auto &s: mySample)
        --s;
}

// Most of the code for rleCpp was obtained
// from Hadley Wickham's article titled 
// "High Performance functions with Rcpp"
// found: http://adv-r.had.co.nz/Rcpp.html
std::vector<int> rleCpp(const std::vector<int> &x) {
    std::vector<int> lengths;
    int prev = x[0];
    std::size_t i = 0;
    lengths.push_back(1);
    
    for(auto it = x.cbegin() + 1; it != x.cend(); ++it) {
        if (prev == *it) {
            ++lengths[i];
        } else {
            lengths.push_back(1);
            prev = *it;
            ++i;
        }
    }
    
    return lengths;
}

double NumPermsWithRep(const std::vector<int> &v) {
    std::vector<int> myLens = rleCpp(v);
    std::sort(myLens.begin(), myLens.end(), std::greater<int>());
    
    const int myMax = myLens[0];
    const int numUni = myLens.size();
    double result = 1;
    
    for (int i = v.size(); i > myMax; --i)
        result *= i;
    
    if (numUni > 1)
        for (int i = 1; i < numUni; ++i)
            for (int j = 2; j <= myLens[i]; ++j)
                result /= j;
    
    return result;
}

double NumPermsNoRep(int n, int k) {
    const double dblN = static_cast<double>(n);
    const double nMinusK = dblN - static_cast<double>(k);
    double result = 1.0;
    
    for (double i = n; i > nMinusK; --i)
        result *= i;
    
    return result;
}

// Returns number of k-combinations from n elements.
// Mathematically speaking, we have: n!/(k!*(n-k)!)
double nChooseK(int n, int k) {
    
    if (k == n || k == 0)
        return 1.0;
    
    double nCk = 1;
    
    for (double i = (n - k + 1), d = 1; d <= k; ++i, ++d) {
        nCk *= i;
        nCk /= d;
    }
    
    return std::round(nCk);
}

double NumCombsWithRep(int n, int r) {
    return nChooseK(n + r - 1, r);
}

// The resulting vector, "triangleVec" resembles triangle numbers. In
// fact, this vector is obtained in a very similar method as generating
// triangle numbers, albeit in a repeating fashion. Two things to keep
// in mind is that we can't guarantee the following:
//      1) the repetition of each element is greater than or equal to n
//      2) that the repetition of each element isn't the same

double MultisetCombRowNumFast(int n, int r, const std::vector<int> &Reps) {
    
    if (r < 1 || n <= 1)
        return 1.0;
    
    if (r == n)
        if (std::accumulate(Reps.cbegin(), Reps.cend(), 0) == n)
            return 1.0;
        
    const int r1 = r + 1;
    std::vector<double> triangleVec(r1);
    std::vector<double> temp(r1);
    
    int myMax = r1;
    
    if (myMax > Reps[0] + 1)
        myMax = Reps[0] + 1;
    
    for (int i = 0; i < myMax; ++i)
        triangleVec[i] = temp[i] = 1;
    
    --myMax;
    int ind = 1;
    
    for (; myMax < r; ++ind) {
        int myMin = std::min(Reps[ind], r);
        
        for (int i = 1; i <= myMin; ++i)
            triangleVec[i] += triangleVec[i - 1];
        
        myMin = std::min(Reps[ind] + myMax, r);
        int j = 0;
        
        for (int i = (Reps[ind] + 1); i <= myMin; ++i, ++j) {
            triangleVec[i] += triangleVec[i - 1];
            triangleVec[i] -= temp[j];
            temp[j] = triangleVec[j];
        }
        
        for (; j <= myMin; ++j)
            temp[j] = triangleVec[j];
        
        myMax = myMin;
    }
    
    const int n1 = n - 1;
    
    for (; ind < n1; ++ind) {
        double t = triangleVec[r];
        const int s = std::min(Reps[ind], r);
        
        for (int i = 1; i <= s; ++i)
            triangleVec[r] += triangleVec[r - i];
        
        double mySum = triangleVec[r];
        
        for (int i = r - 1; i >= s; --i) {
            mySum -= t;
            t = triangleVec[i];
            mySum += triangleVec[i - s];
            triangleVec[i] = mySum;
        }
        
        for (int i = s - 1; i > 0; --i) {
            mySum -= t;
            t = triangleVec[i];
            triangleVec[i] = mySum;
        }
    }
    
    if (ind < n) {
        const int myMin2 = std::min(Reps[n1], r);
        
        for (int i = 1; i <= myMin2; ++i)
            triangleVec[r] += triangleVec[r - i];
    }
    
    return triangleVec[r];
}

// The algorithm below is credited to Randy Lai,
// author of arrangements and iterpc. It is much
// faster than the original naive approach whereby
// we create all combinations of the multiset, then
// subsequently count the number of permutations
// of each of those combinations.
double MultisetPermRowNum(int n, int r, const std::vector<int> &Reps) {
    
    if (n < 2 || r < 1)
        return 1.0;
    
    int sumFreqs = std::accumulate(Reps.cbegin(), Reps.cend(), 0);
    
    if (r > sumFreqs)
        return 0.0;
    
    const int n1 = n - 1;
    const int maxFreq = *std::max_element(Reps.cbegin(), Reps.cend());
    const int myMax = (r < maxFreq) ? (r + 1) : (maxFreq + 1);
    
    // factorial(171)
    // [1] Inf 
    // factorial(170)
    // [1] 7.257416e+306
    if (myMax > 170 || r > 170) {
        mpz_t result;
        mpz_init(result);
        MultisetPermRowNumGmp(result, n, r, Reps);
        
        if (mpz_cmp_d(result, Significand53) > 0)
            return std::numeric_limits<double>::infinity();
        else
            return mpz_get_d(result);
    }
    
    std::vector<int> seqR(r);
    std::iota(seqR.begin(), seqR.end(), 1);
    const double prodR = std::accumulate(seqR.cbegin(), seqR.cend(), 
                                         1.0, std::multiplies<double>());
    
    // Create seqeunce from 1 to myMax, then add another
    // 1 at the front... equivalent to c(1, 1:myMax)
    std::vector<double> cumProd(myMax), resV(r + 1, 0.0);
    std::iota(cumProd.begin(), cumProd.end(), 1);
    
    cumProd.insert(cumProd.begin(), 1);
    std::partial_sum(cumProd.begin(), cumProd.end(), 
                     cumProd.begin(), std::multiplies<double>());
    
    double numPerms = 0.0;
    int myMin = std::min(r, Reps[0]);
    
    for (int i = 0; i <= myMin; ++i)
        resV[i] = prodR / cumProd[i];
    
    for (int i = 1; i < n1; ++i) {
        for (int j = r; j > 0; --j) {
            myMin = std::min(j, Reps[i]);
            numPerms = 0;
            
            for (int k = 0; k <= myMin; ++k)
                numPerms += resV[j - k] / cumProd[k];
            
            resV[j] = numPerms;
        }
    }
    
    myMin = std::min(r, Reps[n1]);
    numPerms = 0;
    
    for (int i = 0; i <= myMin; ++i)
        numPerms += resV[r - i] / cumProd[i];
    
    return numPerms;
}

// This function will be used in the main function to
// determine whether gmp analogs are needed as the fast
// algorithm above could potentially produce negative
// results because of issues with double precision
double MultisetCombRowNum(int n, int r, const std::vector<int> &Reps) {
    
    if (r < 1 || n <= 1)
        return 1;
    
    int i, k, j, myMax, r1 = r + 1;
    std::vector<double> triangleVec(r1);
    std::vector<double> temp(r1);
    
    myMax = r1;
    
    if (myMax > Reps[0] + 1)
        myMax = Reps[0] + 1;
    
    for (i = 0; i < myMax; ++i)
        triangleVec[i] = 1;
    
    temp = triangleVec;
    
    for (k = 1; k < n; ++k) {
        for (i = r; i > 0; --i) {
            myMax = i - Reps[k];
            
            if (myMax < 0) {myMax = 0;}
            
            double tempSum = 0;
            
            for (j = myMax; j <= i; ++j)
                tempSum += triangleVec[j];
            
            temp[i] = tempSum;
        }
        
        triangleVec = temp;
    }
    
    return triangleVec[r];
}

double GetComputedRows(bool IsMultiset, bool IsComb, bool IsRep, int n, int &m, SEXP Rm,
                       int lenFreqs, std::vector<int> &freqs, std::vector<int> &Reps) {
    
    double computedRows = 0;
    
    if (IsMultiset) {
        if (n != lenFreqs)
            Rcpp::stop("the length of freqs must equal the length of v");
        
        if (m > static_cast<int>(freqs.size()))
            m = freqs.size();
        
        if (IsComb) {
            computedRows = MultisetCombRowNum(n, m, Reps);
        } else {
            if (Rf_isNull(Rm) || m == static_cast<int>(freqs.size()))
                computedRows = NumPermsWithRep(freqs);
            else
                computedRows = MultisetPermRowNum(n, m, Reps);
        }
    } else {
        if (IsRep) {
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
    
    return computedRows;
}
