#include <Rcpp.h>
#include <cmath>
#include <RcppThread.h>
#include "PrimesPolRho.h"
#include "CleanConvert.h"

template <typename typeReturn>
void getPrimeFactors(int64_t &t, std::vector<typeReturn> &factors);

const double my63Max = std::pow(2, 63);
const int64_t Sqrt63Max = static_cast<int64_t>(std::sqrt(static_cast<double>(my63Max)));

/* Number of Miller-Rabin tests to run when not proving primality. */
constexpr int MR_REPS = 25;

constexpr int64_t FirstOmittedPrime = 3989;
constexpr std::size_t pDiffSize = sizeof(primesDiffPR) / sizeof(primesDiffPR[0]);

template <typename typeReturn>
void FactorTrialDivision(int64_t &t, std::vector<typeReturn> &factors) {
    
    int p = 3;
    
    while ((t & 1) == 0) {
        factors.push_back(2);
        t >>= 1;
    }
    
    for (std::size_t i = 1; i < pDiffSize;) {
        if ((t % p) != 0) {
            p += primesDiffPR[i++];
            if (t < (p * p))
                break;
        } else {
            t /= p;
            factors.push_back(static_cast<typeReturn>(p));
        }
    }
}

int64_t PositiveMod(int64_t i, int64_t n) {
    return ((i % n) + n) % n;
}

int64_t myGCD(int64_t u, int64_t v) {
    
    u = PositiveMod(u, v);
    
    while (v != 0) {
        int64_t r = u % v;
        u = v;
        v = r;
    }
    
    return u;
}

int64_t ProdBigMod(int64_t x1_i64, int64_t x2_i64, int64_t p_i64) {
    
    double prodX = static_cast<double>(x1_i64) * static_cast<double>(x2_i64);
    int64_t result = 0;
    
    if (prodX < static_cast<double>(p_i64)) {
        result = prodX;
    } else if (p_i64 < Sqrt63Max || prodX < my63Max) {
        result = (x1_i64 * x2_i64) % p_i64;
    } else {
        int64_t nChunks = 1;
        int64_t chunk;
        double part1 = my63Max;
        
        while (part1 >= my63Max) {
            int64_t cSize = static_cast<int64_t>(my63Max / x1_i64);
            chunk = (x1_i64 * cSize) % p_i64;
            nChunks = x2_i64 / cSize;
            int64_t part2 = ((x2_i64 - cSize * nChunks) * x1_i64) % p_i64;
            part1 = static_cast<double>(nChunks) * static_cast<double>(chunk);
            x1_i64 = chunk;
            x2_i64 = nChunks;
            result = (result + part2) % p_i64;
        }
        
        int64_t part1_i64 = (nChunks * chunk) % p_i64;
        result = (part1_i64 + result) % p_i64;
    }
    
    return result;
}

int64_t ExpBySquaring(int64_t x, int64_t n, int64_t p) {
    
    int64_t result;
    
    if (n == 1) {
        result = PositiveMod(x, p);
    } else if (n % 2 == 0) {
        result = ExpBySquaring(ProdBigMod(x, x, p), n / 2, p);
    } else {
        result = ProdBigMod(x, ExpBySquaring(ProdBigMod(x, x, p), (n - 1) / 2, p), p);
    }
    
    return result;
}

int MillerRabin(int64_t n, int64_t nm1, int64_t x,
                int64_t& y, int64_t q, uint64_t k) {
    
    y = ExpBySquaring(x, q, n);
    if (y == 1 || y == nm1)
        return 1;
    
    for (std::size_t i = 1; i < k; ++i) {
        y = ExpBySquaring(y, 2, n);
        if (y == nm1)
            return 1;
        if (y == 1)
            return 0;
    }
    
    return 0;
}

int IsPrime(int64_t n) {
    int k, primeTestReturn;
    int64_t nm1, q, tmp, a;
    
    std::vector<int64_t> factors;
    
    if (n < 2)
        return 0;
    
    /* We have already casted out small primes. */
    if (n < FirstOmittedPrime * FirstOmittedPrime)
        return 1;
    
    /* Precomputation for Miller-Rabin.  */
    q = nm1 = n - 1;
    
    /* Find q and k, where q is odd and n = 1 + 2**k * q.  */
    k = 0;
    while ((q & 1) == 0) {
        q /= 2;
        ++k;
    }
    
    a = 2;
    
    /* Perform a Miller-Rabin test, finds most composites quickly.  */
    if (!MillerRabin (n, nm1, a, tmp, q, (uint64_t) k)) {
        primeTestReturn = 0;
        goto ret2;
    }
    
    /* Factor n-1 for Lucas.  */
    tmp = nm1;
    getPrimeFactors(tmp, factors);
    
    /* Loop until Lucas proves our number prime, or Miller-Rabin proves our
    number composite.  */
    for (std::size_t r = 0; r < pDiffSize; ++r) {
        
        primeTestReturn = 1;
        
        for (std::size_t i = 0; i < factors.size() && primeTestReturn; ++i) {
            tmp = nm1 / factors[i];
            tmp = ExpBySquaring(a, tmp, n);
            primeTestReturn = (tmp != 1);
        }
        
        if (primeTestReturn)
            goto ret1;
        
        a += primesDiffPR[r];	/* Establish new base.  */
            
        if (!MillerRabin(n, nm1, a, tmp, q, (uint64_t) k)) {
            primeTestReturn = 0;
            goto ret1;
        }
    }
    
    Rcpp::stop("Lucas prime test failure. This should not happen");
    ret1:
        factors.resize(0);
        
    ret2:
        return primeTestReturn;
}

void PollardRhoMpzT(mpz_t n, std::size_t a, std::vector<double> &factors) {
    
    mpz_t x, z, y, P, t;
    std::size_t k, q;
    
    mpz_init(t);
    mpz_init_set_si(y, 2);
    mpz_init_set_si(x, 2);
    mpz_init_set_si(z, 2);
    mpz_init_set_ui(P, 1);
    k = q = 1;
    
    while(mpz_cmp_ui(n, 1) != 0) {
        for (;;) {
            do {
                mpz_mul(t, x, x);
                mpz_mod(x, t, n);
                mpz_add_ui(x, x, a);
                
                mpz_sub(t, z, x);
                mpz_mul(t, P, t);
                mpz_mod(P, t, n);
                
                if (k % 32 == 1) {
                    mpz_gcd(t, P, n);
                    if (mpz_cmp_ui(t, 1) != 0)
                        goto factor_found;
                    mpz_set(y, x);
                }
                
            } while (--k != 0);
            
            mpz_set(z, x);
            k = q;
            q <<= 1;
            
            for (std::size_t i = 0; i < k; ++i) {
                mpz_mul(t, x, x);
                mpz_mod(x, t, n);
                mpz_add_ui(x, x, a);
            }
            
            mpz_set(y, x);
        }
        
        factor_found:
            do {
                mpz_mul(t, y, y);
                mpz_mod(y, t, n);
                mpz_add_ui(y, y, a);
                
                mpz_sub(t, z, y);
                mpz_gcd(t, t, n);
                
            } while (mpz_cmp_ui(t, 1) == 0);
        
        mpz_divexact(n, n, t);	/* divide by t, before t is overwritten */

        if (mpz_probab_prime_p(t, MR_REPS) == 0) {
            PollardRhoMpzT(t, a + 1, factors);
        } else {
            double dblT = mpz_get_d(t);
            factors.push_back(dblT);
            
            while (mpz_divisible_p(n, t)) {
                mpz_divexact(n, n, t);
                factors.push_back(dblT);
            }
        }
        
        if (mpz_probab_prime_p(n, MR_REPS) != 0) {
            factors.push_back(mpz_get_d(n));
            break;
        }
        
        mpz_mod(x, x, n);
        mpz_mod(z, z, n);
        mpz_mod(y, y, n);
    }
    
    mpz_clear(P);
    mpz_clear(t);
    mpz_clear(z);
    mpz_clear(x);
    mpz_clear(y);
}

void PollardRho(int64_t n, int64_t a, std::vector<int>& factors) {
    
    int64_t x, z, y, P, t;
    std::size_t k = 1u;
    std::size_t q = 1u;
    
    y = x = z = 2;
    P = 1;
    
    while (n != 1) {
        for (;;) {
            do {
                x *= x;
                x %= n;
                x += a;
                
                t = z - x;
                t = PositiveMod(t, n);
                P *= t;
                P %= n;
                
                if (k % 32 == 1) {
                    t = myGCD(P, n);
                    if (t != 1)
                        goto factor_found;
                    y = x;
                }
                
            } while (--k != 0);
            
            z = x;
            k = q;
            q <<= 1;
            
            for (std::size_t i = 0; i < k; ++i) {
                x *= x;
                x %= n;
                x += a;
            }
            
            y = x;
        }
        
        factor_found:
            do {
                y *= y;
                y %= n;
                y += a;
                t = myGCD(z - y, n);
                
            } while (t == 1);
        
        n /= t;	/* divide by t, before t is overwritten */
        
        if (!IsPrime(t)) {
            PollardRho(t, a + 1, factors);
        } else {
            factors.push_back(static_cast<int>(t));
            
            while ((n % t) == 0) {
                n /= t;
                factors.push_back(static_cast<int>(t));
            }
        }
        
        if (IsPrime(n)) {
            factors.push_back(static_cast<int>(n));
            break;
        }
        
        x %= n;
        z %= n;
        y %= n;
    }
}

template <typename typeReturn>
void getPrimeFactors(int64_t& t, std::vector<typeReturn>& factors) {
    FactorTrialDivision(t, factors);
    
    if (t > 1) {
        if (t < std::numeric_limits<int>::max()) {
            if (IsPrime(t)) {
                factors.push_back(t);
            } else {
                std::vector<int> intFactors;
                PollardRho(t, 1, intFactors);
                factors.insert(factors.end(), intFactors.cbegin(), intFactors.cend());
            }
        } else {
            mpz_t bigT;
            mpz_init(bigT);
            mpz_set_d(bigT, static_cast<double>(t));
            
            if (mpz_probab_prime_p(bigT, MR_REPS)) {
                factors.push_back(t);
            } else {
                std::vector<double> dblFactors;
                PollardRhoMpzT(bigT, 1, dblFactors);
                factors.insert(factors.end(),
                               std::make_move_iterator(dblFactors.cbegin()),
                               std::make_move_iterator(dblFactors.cend()));
            }
            
            mpz_clear(bigT);
        }
    }
    
    std::sort(factors.begin(), factors.end());
}

template <typename typeReturn>
void PrimeFacList(std::size_t m, std::size_t n, std::vector<double> myNums,
                  std::vector<std::vector<typeReturn>> &MyPrimeList) {
    
    for (std::size_t i = m; i < n; ++i) {
        std::vector<typeReturn> factors;
        
        int64_t mPass = static_cast<int64_t>(myNums[i]);
        
        if (mPass < 0) {
            mPass = std::abs(mPass);
            factors.push_back(-1);
        }
        
        if (mPass > 0) {
            getPrimeFactors(mPass, factors);
            MyPrimeList[i] = factors;
        }
    }
}

template <typename typeReturn>
std::vector<typeReturn> Factorize(std::vector<typeReturn> &factors) {
    
    std::size_t n = factors.size();
    
    if (n == 1) {
        std::vector<typeReturn> primeReturn(2, 1);
        primeReturn[1] = factors[0];
        return primeReturn;
    } else {
        std::vector<std::size_t> lengths;
        typeReturn prev = factors[0];
        
        std::size_t numUni = 0;
        std::vector<typeReturn> uniFacs(n);
        uniFacs[0] = factors[0];
        lengths.push_back(1);
        
        for(auto it = factors.cbegin() + 1; it < factors.cend(); ++it) {
            if (prev == *it) {
                ++lengths[numUni];
            } else {
                ++numUni;
                prev = *it;
                lengths.push_back(1);
                uniFacs[numUni] = *it;
            }
        }
        
        std::size_t numFacs = 1;
        
        for (std::size_t i = 0; i <= numUni; ++i)
            numFacs *= (lengths[i] + 1);
        
        std::vector<typeReturn> myFacs(numFacs);

        for (std::size_t i = 0; i <= lengths[0]; ++i)
            myFacs[i] = static_cast<typeReturn>(std::pow(uniFacs[0], i));
        
        if (numUni > 0) {
            std::size_t fSz = 1;
            typeReturn temp;
            
            for (std::size_t j = 1; j <= numUni; ++j) {
                fSz *= (lengths[j - 1] + 1);
                for (std::size_t i = 1; i <= lengths[j]; ++i) {
                    for (std::size_t k = 0, ind = (i * fSz); k < fSz; ++k, ++ind) {
                        temp = static_cast<typeReturn>(std::pow(uniFacs[j], i));
                        temp *= myFacs[k];
                        myFacs[ind] = temp;
                    }
                }
            }
        }
        
        std::sort(myFacs.begin(), myFacs.end());
        return myFacs;
    }
}

template <typename typeReturn>
void FactorList(std::size_t m, std::size_t n, std::vector<double> &myNums,
                std::vector<std::vector<typeReturn>> &MyDivList) {
    
    bool isNegative = false;
    
    for (std::size_t j = m; j < n; ++j) {
        std::vector<typeReturn> myDivisors;
        int64_t mPass = static_cast<int64_t>(myNums[j]);
        
        if (mPass < 0) {
            mPass = std::abs(mPass);
            isNegative = true;
        } else {
            isNegative = false;
        }
        
        if (mPass > 1) {
            std::vector<typeReturn> factors;
            getPrimeFactors(mPass, factors);
            myDivisors = Factorize<typeReturn>(factors);
            
            if (isNegative) {
                const std::size_t facSize = myDivisors.size();
                std::vector<typeReturn> tempInt(2 * facSize);
                std::size_t posInd = facSize, negInd = facSize - 1;
                
                for (std::size_t i = 0; i < facSize; ++i, ++posInd, --negInd) {
                    tempInt[negInd] = -1 * myDivisors[i];
                    tempInt[posInd] = myDivisors[i];
                }
                
                myDivisors = tempInt;
            }
        } else {
            if (isNegative)
                myDivisors.push_back(-1);
            if (mPass > 0)
                myDivisors.push_back(1);
        }
        
        MyDivList[j] = myDivisors;
    }
}

void IsPrimeVec(std::size_t m, std::size_t n, std::vector<double> &myNums,
                Rcpp::LogicalVector &primeTest) {
    
    mpz_t testMpzt;
    mpz_init(testMpzt);
    
    for (std::size_t j = m; j < n; ++j) {
        int64_t testVal = static_cast<int64_t>(myNums[j]);
        
        if (testVal == 1) {
            primeTest[j] = false;
        } else if ((testVal & 1) == 0) {
            if (testVal > 2)
                primeTest[j] = false;
        } else {
            int p = 3;
            for (std::size_t i = 1; i < pDiffSize;) {
                if ((testVal % p) != 0) {
                    p += primesDiffPR[i++];
                    if (testVal < (p * p))
                        break;
                } else {
                    if (testVal > p)
                        primeTest[j] = false;
                    break;
                }
            }
        }
        
        if (primeTest[j]) {
            // 1e9 was determined empirically. Creating an mpz_t and calling
            // mpz_probab_prime_p isn't as efficient as calling a similar
            // function that only deals with primitive types.
            if (myNums[j] < 1000000000) {
                primeTest[j] = IsPrime(testVal);
            } else {
                mpz_set_d(testMpzt, myNums[j]);
                
                if (mpz_probab_prime_p(testMpzt, MR_REPS) == 0)
                    primeTest[j] = false;
            }
        }
    }
    
    mpz_clear(testMpzt);
}

template <typename typeReturn>
void PollardRhoMaster(std::vector<double> &myNums, typeReturn myMax, bool bPrimeFacs,
                      bool bAllFacs, std::vector<std::vector<typeReturn>> &MyList,
                      Rcpp::LogicalVector &primeTest, std::size_t myRange, 
                      int nThreads = 1, int maxThreads = 1) {
    
    bool Parallel = false;
    std::size_t m = 0u;
    
    if (nThreads > 1 && myRange > 1 && maxThreads > 1) {
        Parallel = true;
        if (nThreads > maxThreads) {nThreads = maxThreads;}
        if ((myRange / nThreads) < 1) {nThreads = myRange;}
    }
    
    if (Parallel) {
        RcppThread::ThreadPool pool(nThreads);
        const std::size_t chunkSize = myRange / nThreads;
        std::size_t n = chunkSize - 1;
        
        for (int j = 0; j < (nThreads - 1); m = n, n += chunkSize, ++j) {
            if (bPrimeFacs)
                pool.push(std::cref(PrimeFacList<typeReturn>), m, n, std::ref(myNums), std::ref(MyList));
            else if (bAllFacs)
                pool.push(std::cref(FactorList<typeReturn>), m, n, std::ref(myNums), std::ref(MyList));
            else
                pool.push(std::cref(IsPrimeVec), m, n, std::ref(myNums), std::ref(primeTest));
        }
        
        if (bPrimeFacs)
            pool.push(std::cref(PrimeFacList<typeReturn>), m, myRange, std::ref(myNums), std::ref(MyList));
        else if (bAllFacs)
            pool.push(std::cref(FactorList<typeReturn>), m, myRange, std::ref(myNums), std::ref(MyList));
        else
            pool.push(std::cref(IsPrimeVec), m, myRange, std::ref(myNums), std::ref(primeTest));
        
        pool.join();
        
    } else {
        if (bPrimeFacs)
            PrimeFacList(m, myRange, myNums, MyList);
        else if (bAllFacs)
            FactorList(m, myRange, myNums, MyList);
        else
            IsPrimeVec(m, myRange, myNums, primeTest);
    }
}

template <typename typeReturn>
SEXP TheGlue(std::vector<double> &myNums, typeReturn myMax, bool bPrimeFacs,
             bool bAllFacs, bool keepNames, int nThreads, int maxThreads) {
    
    std::size_t myRange = myNums.size();
    Rcpp::LogicalVector tempVec;
    
    if (bPrimeFacs) {
        if (myRange == 1) {
            int64_t mPass = static_cast<int64_t>(myNums[0]);
            if (mPass == 0) {return Rcpp::IntegerVector();}
            std::vector<typeReturn> factors;
            
            if (mPass < 0) {
                mPass = std::abs(mPass);
                factors.push_back(-1);
            }
            
            getPrimeFactors(mPass, factors);
            return Rcpp::wrap(factors);
        } else {
            std::vector<std::vector<typeReturn>> 
                MyPrimeList(myRange, std::vector<typeReturn>());
            
            PollardRhoMaster(myNums, myMax, bPrimeFacs, bAllFacs,
                             MyPrimeList, tempVec, myRange, nThreads, maxThreads);
            
            Rcpp::List myList = Rcpp::wrap(MyPrimeList);
            if (keepNames)
                myList.attr("names") = myNums;
            
            return myList;
        }
    } else if (bAllFacs) {
        if (myRange == 1) {
            int64_t mPass = static_cast<int64_t>(myNums[0]);
            std::vector<typeReturn> factors;
            std::vector<typeReturn> myDivisors;
            bool isNegative = false;
            
            if (mPass < 0) {
                mPass = std::abs(mPass);
                isNegative = true;
            }
            
            if (mPass > 1) {
                getPrimeFactors(mPass, factors);
                myDivisors = Factorize<typeReturn>(factors);
                
                if (isNegative) {
                    const std::size_t facSize = myDivisors.size();
                    std::vector<typeReturn> negPosFacs(2 * facSize);
                    std::size_t posInd = facSize, negInd = facSize - 1;
                    
                    for (std::size_t i = 0; i < facSize; ++i, ++posInd, --negInd) {
                        negPosFacs[negInd] = -1 * myDivisors[i];
                        negPosFacs[posInd] = myDivisors[i];
                    }
                    
                    return Rcpp::wrap(negPosFacs);
                }
                
                return Rcpp::wrap(myDivisors);
            } else  {
                if (isNegative)
                    myDivisors.push_back(-1);
                if (mPass > 0)
                    myDivisors.push_back(1);
                
                return Rcpp::wrap(myDivisors);
            }
        }
        
        std::vector<std::vector<typeReturn>> 
            MyPrimeList(myRange, std::vector<typeReturn>());
        
        PollardRhoMaster(myNums, myMax, bPrimeFacs, bAllFacs,
                         MyPrimeList, tempVec, myRange, nThreads, maxThreads);
        
        Rcpp::List myList = Rcpp::wrap(MyPrimeList);
        if (keepNames)
            myList.attr("names") = myNums;
        
        return myList;
    } else {
        Rcpp::LogicalVector isPrimeVec(myRange, true);
        std::vector<std::vector<typeReturn>> tempList;
        
        PollardRhoMaster(myNums, myMax, bPrimeFacs, bAllFacs,
                         tempList, isPrimeVec, myRange, nThreads, maxThreads);
        if (keepNames)
            isPrimeVec.attr("names") = myNums;
        
        return isPrimeVec;
    }
}

// [[Rcpp::export]]
SEXP PollardRhoContainer(SEXP Rv, SEXP RNamed, bool bPrimeFacs,
                         bool bAllFacs, SEXP RNumThreads, int maxThreads) {
    
    std::vector<double> myNums;
    bool isNamed = CleanConvert::convertLogical(RNamed, "namedList");

    if (bPrimeFacs || bAllFacs)  // numOnly = true, checkWhole = true, negPoss = true
        CleanConvert::convertVector(Rv, myNums, "v", true, true, true);
    else
        CleanConvert::convertVector(Rv, myNums, "v");

    double myMax = *std::max_element(myNums.cbegin(), myNums.cend());
    double myMin = *std::min_element(myNums.cbegin(), myNums.cend());
    
    if (std::abs(myMin) > myMax)
        myMax = std::abs(myMin);
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads");
    
    if (myMax > std::numeric_limits<int>::max()) {
        return TheGlue(myNums, myMax, bPrimeFacs, bAllFacs, isNamed, nThreads, maxThreads);
    } else {
        int intMax = static_cast<int>(myMax);
        return TheGlue(myNums, intMax, bPrimeFacs, bAllFacs, isNamed, nThreads, maxThreads);
    }
}

