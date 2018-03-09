#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "PrimesPolRho.h"
#include "PollardRho.h"
#include <stdint.h>
using namespace Rcpp;

/* Prove primality or run probabilistic tests.  */
int FlagProvePrimality = 1;

const double Significand53 = 9007199254740991.0;
const double myMax = std::pow(2, 60);
const double my64Max = std::pow(2, 63);
const int64_t Sqrt64Max = (int64_t) std::sqrt((double) my64Max);

#define pDiffSize (sizeof(primesDiffPR) / sizeof(primesDiffPR[0]))

/* Number of Miller-Rabin tests to run when not proving primality. */
#define MR_REPS 25

template <typename typeReturn>
void FactorTrialDivision (int64_t& t,
                          std::vector<typeReturn>& factors) {
    while ((t & 1) == 0) {
        factors.push_back(2);
        t /= 2;
    }
    
    int p = 3;
    for (std::size_t i = 1; i < pDiffSize;) {
        if ((t % p) != 0) {
            p += primesDiffPR[i++];
            if (t < (p * p))
                break;
        } else {
            t /= p;
            factors.push_back((typeReturn) p);
        }
    }
}

inline int64_t PositiveMod(int64_t i, int64_t n) {
    return (i % n + n) % n;
}

int64_t ProdBigMod (int64_t x1_i64,
                         int64_t x2_i64, int64_t p_i64) {
    x1_i64 = PositiveMod(x1_i64, p_i64);
    x2_i64 = PositiveMod(x2_i64, p_i64);
    
    double prodX = (double) (x1_i64) * (double) (x2_i64);
    int64_t result = 0;
    
    if (prodX < (double) p_i64) {
        result = prodX;
    } else if (p_i64 < Sqrt64Max || prodX < my64Max) {
        result = (x1_i64 * x2_i64) % p_i64;
    } else {
        int64_t part2, numChunkMods = 1;
        int64_t chunkSize, chunkMod;
        double part1 = (double) myMax;
        
        while (part1 >= (double) myMax) {
            chunkSize = myMax / x1_i64;
            chunkMod = (x1_i64 * chunkSize) % p_i64;
            numChunkMods = x2_i64 / chunkSize;
            part2 = ((x2_i64 - chunkSize * numChunkMods) * x1_i64) % p_i64;
            part1 = (double) (numChunkMods) * (double) (chunkMod);
            x1_i64 = chunkMod;
            x2_i64 = numChunkMods;
            result = (result + part2) % p_i64;
        }
        
        int64_t part1_i64 = (numChunkMods * chunkMod) % p_i64;
        result = (part1_i64 + result) % p_i64;
    }
    
    return (double) (result);
}

int64_t ExpBySquaring(int64_t x, int64_t n, int64_t p) {
    int64_t result;
    if (n == 1) {
        result = PositiveMod(x, p);
    } else if (n % 2 == 0) {
        result = ExpBySquaring(ProdBigMod(x, x, p), n/2, p);
    } else {
        result = ProdBigMod(x, 
            ExpBySquaring(ProdBigMod(x, x, p), (n - 1)/2, p), p);
    }

    return result;
}

static int64_t myGCD(int64_t u, int64_t v) {
    int64_t r;
    while (v != 0) {
        r = PositiveMod(u, v);
        u = v;
        v = r;
    }
    return u;
}

static int MillerRabin (int64_t n, int64_t nm1,
                        int64_t x, int64_t& y,
                        int64_t q, uint64_t k)
{
    y = ExpBySquaring(x, q, n);
    if (y == 1 || y == nm1)
        return 1;

    for (std::size_t i = 1; i < k; i++) {
        y = ExpBySquaring(y, 2, n);
        if (y == nm1)
            return 1;
        if (y == 1)
            return 0;
    }
    
    return 0;
}

int IsPrime (int64_t n) {
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
        k++;
    }

    a = 2;

    /* Perform a Miller-Rabin test, finds most composites quickly.  */
    if (!MillerRabin (n, nm1, a, tmp, q, (uint64_t) k)) {
        primeTestReturn = 0;
        goto ret2;
    }

    if (FlagProvePrimality) {
        /* Factor n-1 for Lucas.  */
        tmp = nm1;
        getPrimefactors (tmp, factors);
    }

    /* Loop until Lucas proves our number prime, or Miller-Rabin proves our
    number composite.  */
    for (std::size_t r = 0; r < pDiffSize; r++) {

        if (FlagProvePrimality) {
            primeTestReturn = 1;
            for (std::size_t i = 0; i < factors.size() && primeTestReturn; i++) {
                tmp = nm1 / factors[i];
                tmp = ExpBySquaring(a, tmp, n);
                primeTestReturn = (tmp != 1);
            }
        } else {
            /* After enough Miller-Rabin runs, be content. */
            primeTestReturn = (r == (MR_REPS - 1));
        }

        if (primeTestReturn)
            goto ret1;

        a += primesDiffPR[r];	/* Establish new base.  */

        if (!MillerRabin(n, nm1, a, tmp, q, (uint64_t) k)) {
            primeTestReturn = 0;
            goto ret1;
        }
    }

    stop("Lucas prime test failure. This should not happen");
    ret1:
        if (FlagProvePrimality)
            factors.resize(0);

    ret2:
        return primeTestReturn;
}

template <typename typeReturn>
void PollardRho (int64_t n, int64_t a, 
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
            for (i = 0; i < k; i++) {
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
        
        n = n/t;	/* divide by t, before t is overwritten */
        
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
void getPrimefactors (int64_t& t, std::vector<typeReturn>& factors) {
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
List PrimeFacList (std::vector<double> myNums, bool namedList) {
    
    int64_t mPass;
    unsigned int myLen = myNums.size();
    
    std::vector<std::vector<typeReturn> > 
        MyPrimeList(myLen, std::vector<typeReturn>());
    
    for (std::size_t i = 0; i < myLen; i++) {
        std::vector<typeReturn> factors;
        
        mPass = (int64_t) myNums[i];
        
        if (mPass < 0) {
            mPass = std::abs(mPass);
            factors.push_back(-1);
        }
        
        if (mPass > 0) {
            getPrimefactors(mPass, factors);
            MyPrimeList[i] = factors;
        }
    }
    
    Rcpp::List myList = wrap(MyPrimeList);
    if (namedList)
        myList.attr("names") = myNums;
    return myList;
}

// [[Rcpp::export]]
SEXP PrimeFactorsContainer (SEXP Rv, SEXP RNamed) {
    std::vector<double> myNums;
    bool isNamed = as<bool>(RNamed);
    
    switch(TYPEOF(Rv)) {
        case REALSXP: {
            myNums = as<std::vector<double> >(Rv);
            break;
        }
        case INTSXP: {
            myNums = as<std::vector<double> >(Rv);
            break;
        }
        default: {
            stop("v must be of type numeric or integer");
        }
    }
    
    if (myNums.size() > 1) {
        double myMax = *std::max_element(myNums.begin(), myNums.end());
        double myMin = *std::min_element(myNums.begin(), myNums.end());
        
        if (std::abs(myMin) > myMax)
            myMax = std::abs(myMin);
        
        if (myMax > Significand53)
            stop("the abs value of each element must be less than 2^53");
        
        if (myMax > INT_MAX)
            return PrimeFacList<double>(myNums, isNamed);
        else
            return PrimeFacList<int>(myNums, isNamed);
    } else {
        int64_t mPass = myNums[0];
        bool isNegative = false;
        
        if (mPass < 0) {
            mPass = std::abs(mPass);
            isNegative = true;
        }
        
        if (mPass == 0)
            return IntegerVector();
        
        if (mPass > Significand53)
            stop("the abs value of each element must be less than 2^53");
        
        if (mPass > INT_MAX) {
            std::vector<double> factors;
            if (isNegative)
                factors.push_back(-1);
            getPrimefactors(mPass, factors);
            return wrap(factors);
        } else {
            std::vector<int> factors;
            if (isNegative)
                factors.push_back(-1);
            getPrimefactors(mPass, factors);
            return wrap(factors);
        }
    }
}

// [[Rcpp::export]]
SEXP IsPrimeContainer (SEXP Rv, SEXP RNamed) {
    std::vector<double> myNums;
    int64_t testVal;
    bool isNamed = false;
    
    switch(TYPEOF(Rv)) {
        case REALSXP: {
            myNums = as<std::vector<double> >(Rv);
            break;
        }
        case INTSXP: {
            myNums = as<std::vector<double> >(Rv);
            break;
        }
        default: {
            stop("v must be of type numeric or integer");
        }
    }
    
    isNamed = as<bool>(RNamed);
    unsigned int myLen = myNums.size();
    std::vector<bool> primeTest(myLen, true);
    double isWhole;
    
    for (std::size_t j = 0; j < myLen; j++) {
        
        if (myNums[j] <= 0)
            stop("each element must be positive");
        if (myNums[j] > Significand53)
            stop("each element must be less than 2^53");
        if (std::modf(myNums[j], &isWhole) != 0.0)
            primeTest[j] = false;
        
        testVal = (int64_t) myNums[j];
        
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
        
        if (primeTest[j] == 1)
            primeTest[j] = (bool) IsPrime(testVal);
    }
    
    LogicalVector isPrimeVec = wrap(primeTest);
    
    if (isNamed)
        isPrimeVec.attr("names") = myNums;
    
    return isPrimeVec;
}
