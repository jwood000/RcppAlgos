#include <Rcpp.h>
#include <math.h>
#include <algorithm>
#include "PrimesPolRho.h"
#include "PollardRho.h"
using namespace Rcpp;

/* Prove primality or run probabilistic tests.  */
int FlagProvePrimality = 0;

const double Significand53 = 9007199254740991;
const double SqrtSig53 = std::floor(std::sqrt(Significand53));
#define pDiffSize (sizeof(primesDiffPR) / sizeof(primesDiffPR[0]))

/* Number of Miller-Rabin tests to run when not proving primality. */
#define MR_REPS 25

void FactorTrialDivision (double& t,
                          std::vector<double>& factors) {
    unsigned long int p;
    int i = 1;
    
    p = 2;
    for (i = 0; i < pDiffSize;) {
        if (std::fmod(t, p) != 0) {
            p += primesDiffPR[i++];
            if (t < p * p)
                break;
        } else {
            t /= p;
            t = round(t);
            factors.push_back(p);
        }
    }
}

double PositiveMod(double x, double m) {
    if (x < 0)
        x = x + ceil(std::abs(x) / m) * m;
    else if (x > m)
        x = std::fmod(x, m);
    return x;
}

double ProdBigMod(double x1, double x2, double p) {
    double result = 0, prodX;
    
    x1 = PositiveMod(x1, p);
    x2 = PositiveMod(x2, p);
    prodX = x1 * x2;
    
    if (prodX < p) {
        result = prodX;
    } else if (p < SqrtSig53) {
        result = std::fmod(prodX, p);
    } else {
        double numChunkMods, part1 = Significand53;
        double initialChunk, chunkMod, part2;
        
        while (part1 >= Significand53) {
            initialChunk = floor(Significand53 / x1);
            chunkMod = std::fmod(x1 * initialChunk, p);
            numChunkMods = floor(x2 / initialChunk);
            part2 = std::fmod((x2 - initialChunk * numChunkMods) * x1, p);
            part1 = numChunkMods * chunkMod;
            x1 = chunkMod; x2 = numChunkMods;
            result = std::fmod(result + part2, p);
        }
        
        result = std::fmod(part1 + result, p);
    }
    
    return result;
}

double ExpBySquaring(double x, double n, double p) {
    double result;
    if (n == 1) {
        result = PositiveMod(x, p);
    } else if (std::fmod(n, 2) == 0) {
        result = ExpBySquaring(ProdBigMod(x, x, p), n/2, p);
    } else {
        result = ProdBigMod(x, 
            ExpBySquaring(ProdBigMod(x, x, p), (n - 1)/2, p), p);
    }

    return result;
}

static double myGCD(double u, double v) {
    double r;
    while (v != 0) {
        r = PositiveMod(u, v);
        u = v;
        v = r;
    }
    return u;
}

static int MillerRabin (double n, double nm1,
                        unsigned long int x, double& y,
                        double q, unsigned long int k)
{
    unsigned long int i;
    y = ExpBySquaring((double)x, q, n);
    if (y == 1 || y == nm1)
        return 1;

    for (i = 1; i < k; i++) {
        y = ExpBySquaring(y, 2, n);
        if (y == nm1)
            return 1;
        if (y == 1)
            return 0;
    }
    
    return 0;
}

int IsPrime (double n) {
    int k, r, primeTestReturn;
    double nm1, q, tmp;
    unsigned long int a;

    std::vector<double> factors;

    if (n < 2)
        return 0;

    /* We have already casted out small primes. */
    if (n < FirstOmittedPrime * FirstOmittedPrime)
        return 1;

    /* Precomputation for Miller-Rabin.  */
    q = nm1 = n - 1;

    /* Find q and k, where q is odd and n = 1 + 2**k * q.  */
    k = 0;
    while (std::fmod(q, 2) == 0) {
        q /= 2;
        q = round(q);
        k++;
    }

    a = 2;

    /* Perform a Miller-Rabin test, finds most composites quickly.  */
    if (!MillerRabin (n, nm1, a, tmp, q, k)) {
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
    for (r = 0; r < pDiffSize; r++) {
        int i;

        if (FlagProvePrimality) {
            primeTestReturn = 1;
            for (i = 0; i < factors.size() && primeTestReturn; i++) {
                tmp = nm1 / factors[i];
                tmp = ExpBySquaring((double)a, tmp, n);
                primeTestReturn = (tmp != 1);
            }
        } else {
            /* After enough Miller-Rabin runs, be content. */
            primeTestReturn = (r == MR_REPS - 1);
        }

        if (primeTestReturn)
            goto ret1;

        a += primesDiffPR[r];	/* Establish new base.  */

        if (!MillerRabin(n, nm1, a, tmp, q, k)) {
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

void PollardRho (double n, unsigned long a, 
                 std::vector<double>& factors) {
    double x, z, y, P, t;
    unsigned long  k, l, i;
    
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
        
        n = round(n/t);	/* divide by t, before t is overwritten */
        
        if (!IsPrime(t)) {
            PollardRho(t, a + 1, factors);
        } else {
            factors.push_back(t);
        }
        
        if (IsPrime(n)) {
            factors.push_back(n);
            break;
        }
        
        x = PositiveMod(x, n);
        z = PositiveMod(z, n);
        y = PositiveMod(y, n);
    }
}

void getPrimefactors (double& t, std::vector<double>& factors) {
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

// [[Rcpp::export]]
SEXP PrimeFactorsContainer (SEXP n) {
    double m;
    std::vector<double> factors;
    
    switch(TYPEOF(n)) {
        case REALSXP: {
            m = as<double>(n);
            break;
        }
        case INTSXP: {
            m = as<double>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
    if (m < 0) {
        m = std::abs(m);
        factors.push_back(-1);
    }
    
    if (m > 0) {
        getPrimefactors(m, factors);
        return wrap(factors);
    } else {
        IntegerVector trivialReturn;
        return trivialReturn;
    }
}