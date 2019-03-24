#include <Rcpp.h>
#include <cmath>
#include <RcppThread.h>
#include "PrimesPolRho.h"
#include "CleanConvert.h"

// forward declare getPrimeFactors as isPrime needs it and
// PollardRho (which is called by getPrimeFactors) calls isPrime
template <typename typeReturn>
void getPrimeFactors (int64_t &t, std::vector<typeReturn> &factors);

// Prove primality or run probabilistic tests
int FlagProvePrimality = 1;

const double my2Pow62 = std::pow(2, 62);
const double my64Max = std::pow(2, 63);
const int64_t Sqrt64Max = static_cast<int64_t>(std::sqrt(my64Max));

constexpr int64_t FirstOmittedPrime = 4001;
constexpr std::size_t pDiffSize = sizeof(primesDiffPR) / sizeof(primesDiffPR[0]);

/* Number of Miller-Rabin tests to run when not proving primality. */
constexpr std::size_t MR_REPS = 25;

template <typename typeReturn>
void FactorTrialDivision (int64_t& t,
                          std::vector<typeReturn>& factors) {
    while ((t & 1) == 0) {
        factors.push_back(2);
        t >>= 1;
    }
    
    int p = 3;
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

inline int64_t PositiveMod(int64_t i, int64_t n) {
    return ((i % n) + n) % n;
}

int64_t ProdBigMod(int64_t x1_i64, int64_t x2_i64, int64_t p_i64) {
    
    x1_i64 = PositiveMod(x1_i64, p_i64);
    x2_i64 = PositiveMod(x2_i64, p_i64);
    
    double prodX = static_cast<double>(x1_i64) * static_cast<double>(x2_i64);
    int64_t result = 0;
    
    if (prodX < static_cast<double>(p_i64)) {
        result = prodX;
    } else if (p_i64 < (int64_t) Sqrt64Max || prodX < my64Max) {
        result = (x1_i64 * x2_i64) % p_i64;
    } else {
        int64_t part2, numChunkMods = 1;
        int64_t chunkSize, chunkMod;
        double part1 = my2Pow62;
        
        while (part1 >= my2Pow62) {
            chunkSize = static_cast<int64_t>(my2Pow62 / x1_i64);
            chunkMod = (x1_i64 * chunkSize) % p_i64;
            numChunkMods = x2_i64 / chunkSize;
            part2 = ((x2_i64 - (chunkSize * numChunkMods)) * x1_i64) % p_i64;
            part1 = static_cast<double>(numChunkMods) * static_cast<double>(chunkMod);
            x1_i64 = chunkMod;
            x2_i64 = numChunkMods;
            result = (result + part2) % p_i64;
        }
        
        int64_t part1_i64 = (numChunkMods * chunkMod) % p_i64;
        result = (part1_i64 + result) % p_i64;
    }
    
    return static_cast<double>(result);
}

int64_t ExpBySquaring(int64_t x, int64_t n, int64_t p) {
    int64_t result;
    if (n == 1) {
        result = PositiveMod(x, p);
    } else if ((n % 2) == 0) {
        result = ExpBySquaring(ProdBigMod(x, x, p), n / 2, p);
    } else {
        result = ProdBigMod(x, 
            ExpBySquaring(ProdBigMod(x, x, p), (n - 1) / 2, p), p);
    }

    return result;
}

int64_t myGCD(int64_t u, int64_t v) {
    int64_t r;
    while (v != 0) {
        r = PositiveMod(u, v);
        u = v;
        v = r;
    }
    return u;
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
        q >>= 1;
        ++k;
    }

    a = 2;

    /* Perform a Miller-Rabin test, finds most composites quickly.  */
    if (!MillerRabin(n, nm1, a, tmp, q, (uint64_t) k)) {
        primeTestReturn = 0;
        goto ret2;
    }

    if (FlagProvePrimality) {
        /* Factor n-1 for Lucas.  */
        tmp = nm1;
        getPrimeFactors(tmp, factors);
    }

    /* Loop until Lucas proves our number prime, or Miller-Rabin proves our
    number composite.  */
    for (std::size_t r = 0; r < pDiffSize; ++r) {

        if (FlagProvePrimality) {
            primeTestReturn = 1;
            for (std::size_t i = 0; i < factors.size() && primeTestReturn; ++i) {
                tmp = nm1 / factors[i];
                tmp = ExpBySquaring(a, tmp, n);
                primeTestReturn = (tmp != 1);
            }
        } else {
            /* After enough Miller-Rabin runs, be content. */
            primeTestReturn = (r == MR_REPS);
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
        if (FlagProvePrimality)
            factors.resize(0);

    ret2:
        return primeTestReturn;
}

template <typename typeReturn>
void PollardRho(int64_t n, int64_t a, 
                std::vector<typeReturn>& factors) {
    int64_t x, z, y, P, t;
    int64_t  k, l, i;
    
    y = x = z = 2;
    P = k = l = 1;
    
    while (n != 1) {
        for (;;) {
            do {
                x = ProdBigMod(x, x, n);
                x += a;
                t = z - x;
                P = ProdBigMod(P, t, n);
                
                if (k % 32 == 1) {
                    t = myGCD(P, n);
                    if (t != 1)
                        goto factor_found;
                    y = x;
                }
            } while (--k != 0);
            
            z = x;
            k = l;
            l = 2 * l;
            for (i = 0; i < k; ++i) {
                x = ProdBigMod(x, x, n);
                x += a;
            }
            y = x;
        }
        
        factor_found:
            do {
                y = ProdBigMod(y, y, n);
                y += a;
                t = z - y;
                t = myGCD(t, n);
            } while (t == 1);
        
        n /= t;	/* divide by t, before t is overwritten */
        
        if (!IsPrime(t)) {
            PollardRho(t, a + 1, factors);
        } else {
            factors.push_back((typeReturn) t);
        }
        
        if (IsPrime(n)) {
            factors.push_back((typeReturn) n);
            break;
        }
        
        x = PositiveMod(x, n);
        z = PositiveMod(z, n);
        y = PositiveMod(y, n);
    }
}

template <typename typeReturn>
void getPrimeFactors(int64_t& t, std::vector<typeReturn>& factors) {
    FactorTrialDivision(t, factors);
    
    if (t > 1) {
        if (IsPrime(t)) {
            factors.push_back(t);
        } else {
            PollardRho(t, 1, factors);
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
    
    unsigned long int n = factors.size();
    
    if (n == 1) {
        std::vector<typeReturn> primeReturn(2, 1);
        primeReturn[1] = factors[0];
        return primeReturn;
    } else {
        std::vector<unsigned long int> lengths;
        typeReturn prev = factors[0];
        
        unsigned long int numUni = 0;
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
        
        unsigned long int numFacs = 1;
        
        for (std::size_t i = 0; i <= numUni; ++i)
            numFacs *= (lengths[i] + 1);
        
        std::vector<typeReturn> myFacs(numFacs);

        for (std::size_t i = 0; i <= lengths[0]; ++i)
            myFacs[i] = static_cast<typeReturn>(std::pow(uniFacs[0], i));
        
        if (numUni > 0) {
            unsigned long int fSz = 1;
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
    
    int64_t mPass;
    bool isNegative = false;
    
    for (std::size_t j = m; j < n; ++j) {
        std::vector<typeReturn> myDivisors;
        mPass = static_cast<int64_t>(myNums[j]);
        
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
                const unsigned long int facSize = myDivisors.size();
                std::vector<typeReturn> tempInt(2 * facSize);
                unsigned long int posInd = facSize, negInd = facSize - 1;
                
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
        
        if (primeTest[j])
            primeTest[j] = IsPrime(testVal);
    }
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
                    const unsigned long int facSize = myDivisors.size();
                    std::vector<typeReturn> negPosFacs(2 * facSize);
                    unsigned long int posInd = facSize, negInd = facSize - 1;
                    
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

