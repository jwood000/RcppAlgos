#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// Below are two functions for generating prime numbers quickly.
// "AllPrimesC" implements an optimized version of the Sieve of Eratosthenes.
// "RangePrimeC" will return all primes within a given range.

NumericVector AllPrimesCpp (long long int n) {
    std::vector<bool> primes(n+1, true);
    std::vector<long long int> myps;
    int myReserve = (long long int)floor((double)2*n/log((double)n));
    myps.reserve(myReserve);

    long long int lastP = 3;
    long long int fsqr = (long long int)floor((double)sqrt((double)n));
    long long int k, ind, j;

    for (j = 4; j <= n; j += 2) {primes[j] = false;}

    while (lastP <= fsqr) {
        for (j = lastP*lastP; j <= n; j += 2*lastP) {primes[j] = false;}
        k = lastP + 2;
        ind = 2;
        while (!primes[k]) {
            k += 2;
            ind += 2;
        }
        lastP += ind;
    }

    myps.push_back(2);

    for (j = 3; j <= n; j += 2) {
        if (primes[j]) {
            myps.push_back(j);
        }
    }

    return wrap(myps);
}

NumericVector RangePrimesCpp (long long int m, long long int n) {
    long long int limit = n-m;
    std::vector<bool> primes(limit+1, true);
    std::vector<long long int> myps;
    long long int myReserve = (long long int)floor((double)2*limit/log((double)limit));
    myps.reserve(myReserve);

    IntegerVector testPrimes = as<IntegerVector>(AllPrimesCpp((long long int)floor((double)sqrt((double)n))));
    IntegerVector::iterator p, maxP;
    maxP = testPrimes.end();

    long long int j, strt, f_odd, remTest, nMinusM = n- m, powTest;
    double my2 = 2;
    strt = f_odd = m;

    if (m==1) { // Get first odd value
        f_odd=3;
    } else if (f_odd % 2 == 0) {
        f_odd++;
    }

    while (strt % 2 != 0) {strt++;}

    if (m <= 2) {
        strt=4;
        myps.push_back(2);
    }

    for (j = strt; j <= n; j += 2) {primes[j-m] = false;}

    for (p = testPrimes.begin() + 1; p < maxP; p++) {
        powTest = (long long int)pow((double)*p, my2);
        if (m > powTest) {
            remTest = m % *p;
            if (remTest == 0) {
                strt = 0;
            } else {
                strt = *p - remTest;
            }
            if ((m+strt) % 2 == 0) {strt += *p;}
            for (j = strt; j <= nMinusM; j += (*p)*2) {primes[j] = false;}
        } else {
            strt = powTest - m;
            for (j = strt; j <= nMinusM; j += (*p)*2) {primes[j] = false;}
        }
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
    if (mTest < 0) {stop("n must be positive");}
    m = mTest;
    
    std::vector<std::vector<int> > MyPrimeList(m, std::vector<int>());
    std::vector<std::vector<int> >::iterator it2d, itEnd;
    itEnd = MyPrimeList.end();
    IntegerVector primes = as<IntegerVector>(AllPrimesCpp(m));
    IntegerVector::iterator p;
    int i, j, limit, myStep, myMalloc;
    myMalloc = (int)floor((double)log2((double)m));
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
    if (mTest < 0) {stop("n must be positive");}
    m = mTest;

    IntegerVector starterSequence = Rcpp::seq(1, m);
    NumericVector EulerPhis = as<NumericVector>(starterSequence);
    std::vector<int> EulerInt(m);
    IntegerVector primes = as<IntegerVector>(AllPrimesCpp(m));
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
    long long int bound1, bound2, myMax, myMin;
    double b1Test, b2Test;

    switch(TYPEOF(Rb1)) {
        case REALSXP: {
            b1Test = as<double>(Rb1);
            break;
        }
        case INTSXP: {
            b1Test = as<double>(Rb1);
            break;
        }
        default: {
            stop("bound1 must be of type numeric or integer");
        }
    }
    
    if (b1Test < 1 || b1Test > 9007199254740991) {stop("bound1 must be a positive number less than 2^53");}
    bound1 = b1Test;

    if (Rf_isNull(Rb2)) {
        if (bound1 == 1) {
            IntegerVector v;
            return v;
        } else {
            return AllPrimesCpp(bound1);
        }
    } else {
        switch(TYPEOF(Rb2)) {
            case REALSXP: {
                b2Test = as<double>(Rb2);
                break;
            }
            case INTSXP: {
                b2Test = as<double>(Rb2);
                break;
            }
            default: {
                stop("bound2 must be of type numeric or integer");
            }
        }
        if (b2Test < 1 || b2Test > 9007199254740991) {stop("bound2 must be a positive number less than 2^53");}
        bound2 = b2Test;

        if (bound1 > bound2) {
            myMax = bound1;
            myMin = bound2;
        } else {
            myMax = bound2;
            myMin = bound1;
        }

        if (myMax <= 1) {
            IntegerVector z;
            return z;
        }

        if (myMin <= 2) {
            return AllPrimesCpp(myMax);
        } else {
            return RangePrimesCpp(myMin, myMax);
        }
    }
}
