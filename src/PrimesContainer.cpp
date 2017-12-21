#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// Below are two functions for generating prime numbers quickly.
// "AllPrimesC" implements an optimized version of the Sieve of Eratosthenes.
// "RangePrimeC" will return all primes within a given range.

IntegerVector AllPrimesCpp (int n) {
    std::vector<bool> primes(n+1, true);
    std::vector<int> myps;
    int myReserve = (int)floor((double)2*n/log((double)n));
    myps.reserve(myReserve);

    unsigned long int lastP = 3, k, ind, j, uN = n, step;
    unsigned long int fsqr = (unsigned long int)floor(sqrt((double)n));

    for (j = 4; j <= uN; j += 2) {primes[j] = false;}
    myps.push_back(2);
    
    while (lastP <= fsqr) {
        step = 2*lastP;
        for (j = lastP*lastP; j <= uN; j += step) {primes[j] = false;}
        k = lastP + 2;
        ind = 2;
        while (!primes[k]) {
            k += 2;
            ind += 2;
        }
        myps.push_back(lastP);
        lastP += ind;
    }

    for (j = lastP; j <= uN; j += 2) {
        if (primes[j]) {
            myps.push_back(j);
        }
    }

    return wrap(myps);
}

NumericVector RangePrimesCpp (double m, double n) {
    double limit = n-m;
    std::vector<bool> primes(limit+1, true);
    std::vector<double> myps;
    double myReserve = floor(2*limit/log(limit));
    myps.reserve(myReserve);

    IntegerVector testPrimes = AllPrimesCpp((int)floor(sqrt(n)));
    IntegerVector::iterator p, maxP;
    maxP = testPrimes.end();

    double j, strt, f_odd, remTest, nMinusM = n - m, powTest, step;
    strt = f_odd = m;

    if (m==1) { // Get first odd value
        f_odd = 3;
    } else if (fmod(f_odd, 2) == 0.0) {
        f_odd++;
    }

    while (fmod(strt, 2) != 0.0) {strt++;}

    if (m <= 2) {
        strt = 4;
        myps.push_back(2);
    }

    for (j = strt; j <= n; j += 2) {primes[j-m] = false;}

    for (p = testPrimes.begin() + 1; p < maxP; p++) {
        powTest = pow((double)*p, 2);
        step = (*p)*2;
        if (m > powTest) {
            remTest = fmod(m,*p);
            if (remTest == 0.0) {
                strt = 0;
            } else {
                strt = *p - remTest;
            }
            if (fmod(m+strt, 2) == 0.0) {strt += *p;}
        } else {
            strt = powTest - m;
        }
        for (j = strt; j <= nMinusM; j += step) {primes[j] = false;}
    }

    for (j = f_odd; j <= n; j += 2) {
        if (primes[j-m]) {
            myps.push_back(j);
        }
    }

    return wrap(myps);
}

// [[Rcpp::export]]
List PrimeFactorizationListRcpp (SEXP n) {
    int m;
    double mTest;

    switch(TYPEOF(n)) {
        case REALSXP: {
            mTest = as<double>(n);
            break;
        }
        case INTSXP: {
            mTest = as<double>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
    if (mTest > 2147483647) {stop("n must be less than 2^31");}
    
    if (mTest <= 0) {
        stop("n must be positive");
    } else if (mTest < 2) {
        List trivialReturn;
        return trivialReturn;
    }
    
    m = (int)ceil(mTest);
    
    std::vector<std::vector<int> > MyPrimeList(m, std::vector<int>());
    std::vector<std::vector<int> >::iterator it2d, itEnd;
    itEnd = MyPrimeList.end();
    IntegerVector primes = AllPrimesCpp(m);
    IntegerVector::iterator p;
    int i, j, limit, myStep, myMalloc;
    myMalloc = (int)floor(log2((double)m));
    double myLogN = log((double)m);
    
    for (it2d = MyPrimeList.begin(); it2d < itEnd; it2d++) {
        it2d -> reserve(myMalloc);
    }

    for (p = primes.begin(); p < primes.end(); p++) {
        limit = (int)trunc(myLogN/log((double)*p));
        for (i = 1; i <= limit; i++) {
            myStep = (int)pow((double)*p,i);
            for (j = myStep; j <= m; j += myStep) {
                MyPrimeList[j-1].push_back(*p);
            }
        }
    }

    return wrap(MyPrimeList);
}

// [[Rcpp::export]]
IntegerVector EulerPhiSieveRcpp (SEXP n) {
    int m;
    double mTest;

    switch(TYPEOF(n)) {
        case REALSXP: {
            mTest = as<double>(n);
            break;
        }
        case INTSXP: {
            mTest = as<double>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }
    
    if (mTest > 2147483647) {stop("n must be less than 2^31");}
    
    if (mTest <= 0) {
        stop("n must be positive");
    } else if (mTest <= 1) {
        IntegerVector trivialReturn(1, 1);
        return trivialReturn;
    }
    
    m = (int)ceil(mTest);
    
    IntegerVector starterSequence = Rcpp::seq(1, m);
    NumericVector EulerPhis = as<NumericVector>(starterSequence);
    std::vector<int> EulerInt(m);
    IntegerVector primes = AllPrimesCpp(m);
    IntegerVector::iterator p;
    int j;

    for (p = primes.begin(); p < primes.end(); p++) {
        for (j = *p; j <= m; j += *p) {
            EulerPhis[j-1] *= (1.0 - 1.0/(*p));
        }
    }

    for (j = 0; j < m; j++) {EulerInt[j] = round(EulerPhis[j]);}

    return wrap(EulerInt);
}

// [[Rcpp::export]]
SEXP EratosthenesRcpp (SEXP Rb1, SEXP Rb2) {
    double bound1, bound2, myMax, myMin;

    switch(TYPEOF(Rb1)) {
        case REALSXP: {
            bound1 = as<double>(Rb1);
            break;
        }
        case INTSXP: {
            bound1 = as<double>(Rb1);
            break;
        }
        default: {
            stop("bound1 must be of type numeric or integer");
        }
    }

    if (Rf_isNull(Rb2)) {
        if (bound1 <= 0) {stop("n must be positive");}
        
        if (bound1 < 2) {
            IntegerVector v;
            return v;
        } else {
            if (bound1 >= 2147483647) {
                if (bound1 > 2147483647) {
                    stop("n must be less than 2^31");
                } else {
                    // In the AllPrimesCpp function above, the very first
                    // step is to instantiate a logical vector with length
                    // equal to n + 1. This is to avoid unnecessary adding
                    // or subtracting one in order to obtain actual numbers
                    // instead of base zero indices. As a workaround, we
                    // generate all primes up to 2^31 - 3 to avoid exceeding 
                    // the length limit of vectors in C++ and subsequently
                    // append 2^31 - 1, since it is prime.
                    IntegerVector IntegerLimitPrimes = AllPrimesCpp(2147483645);
                    IntegerLimitPrimes.push_back(2147483647);
                    return IntegerLimitPrimes;
                }
            }
            return AllPrimesCpp((int)bound1);
        }
    } else {
        if (bound1 <= 0 || bound1 > 9007199254740991.0) {stop("bound1 must be a positive number less than 2^53");}
        
        switch(TYPEOF(Rb2)) {
            case REALSXP: {
                bound2 = as<double>(Rb2);
                break;
            }
            case INTSXP: {
                bound2 = as<double>(Rb2);
                break;
            }
            default: {
                stop("bound2 must be of type numeric or integer");
            }
        }
        if (bound2 <= 0 || bound2 > 9007199254740991.0) {stop("bound2 must be a positive number less than 2^53");}

        if (bound1 > bound2) {
            myMax = bound1;
            myMin = bound2;
        } else {
            myMax = bound2;
            myMin = bound1;
        }
        
        myMin = ceil(myMin);
        myMax = floor(myMax);

        if (myMax <= 1) {
            IntegerVector z;
            return z;
        }

        if (myMin <= 2) {
            return RangePrimesCpp(1, myMax);
        } else {
            return RangePrimesCpp(myMin, myMax);
        }
    }
}
