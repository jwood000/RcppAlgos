#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// Below are two functions for generating prime numbers quickly.
// "AllPrimesC" implements an optimized version of the Sieve of Eratosthenes.
// "RangePrimeC" will return all primes within a given range.

IntegerVector AllPrimesCpp (int n) {
    std::vector<bool> primes(n+1, true);
    std::vector<int> myps;
    myps.reserve(floor(2*n/log(n)));

    int lastP = 3;
    int fsqr = floor(sqrt(n));
    int k, ind, j;

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

IntegerVector RangePrimesCpp (int m, int n) {
    int limit = n-m;
    std::vector<bool> primes(limit+1, true);
    std::vector<int> myps;
    myps.reserve(floor(2*limit/log(limit)));

    IntegerVector testPrimes = AllPrimesCpp(floor(sqrt(n)));
    IntegerVector::iterator p;

    int j, strt, f_odd;
    strt = f_odd = m;

    if (m==1) { // Get first odd value
        f_odd=3;
    } else if (f_odd % 2 == 0) {
        ++f_odd;
    }

    while (strt % 2 != 0) {++strt;}

    if (m <= 2) {
        strt=4;
        myps.push_back(2);
    }

    for (j = strt; j <= n; j += 2) {primes[j-m] = false;}

    for (p = testPrimes.begin() + 1; p < testPrimes.end(); ++p) {
        if (m > *p) {
            strt = f_odd;
            while (strt % *p != 0) {strt += 2;}
            for (j = strt; j <= n; j += (*p)*2) {primes[j-m] = false;}
        } else {
            strt = f_odd;
            while (strt % *p != 0) {strt += 2;}
            for (j = strt*strt; j <= n; j += (*p)*2) {primes[j-m] = false;}
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

    switch(TYPEOF(n)) {
        case REALSXP: {
            m = as<int>(n);
            break;
        }
        case INTSXP: {
            m = as<int>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }

    std::vector<std::vector<int> > MyPrimeList(m, std::vector<int>());
    IntegerVector primes = AllPrimesCpp(m);
    IntegerVector::iterator p;
    int i, j, limit, myStep;
    double myLogN = log(m);

    for (p = primes.begin(); p < primes.end(); ++p) {
        limit = trunc(myLogN/log(*p));
        for (i = 1; i <= limit; ++i) {
            myStep = pow(*p,i);
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

    switch(TYPEOF(n)) {
        case REALSXP: {
            m = as<int>(n);
            break;
        }
        case INTSXP: {
            m = as<int>(n);
            break;
        }
        default: {
            stop("n must be of type numeric or integer");
        }
    }

    IntegerVector starterSequence = Rcpp::seq(1, m);
    NumericVector EulerPhis = as<NumericVector>(starterSequence);
    std::vector<int> EulerInt(m);
    IntegerVector primes = AllPrimesCpp(m);
    IntegerVector::iterator p;
    int j;

    for (p = primes.begin(); p < primes.end(); ++p) {
        for (j = *p; j <= m; j += *p) {
            EulerPhis[j-1] *= (1.0 - 1.0/(*p));
        }
    }

    for (j = 0; j < m; ++j) {EulerInt[j] = round(EulerPhis[j]);}

    return wrap(EulerInt);
}

// [[Rcpp::export]]
SEXP EratosthenesRcpp (SEXP Rb1, SEXP Rb2) {
    int bound1, bound2, myMax, myMin;

    switch(TYPEOF(Rb1)) {
        case REALSXP: {
            bound1 = as<int>(Rb1);
            break;
        }
        case INTSXP: {
            bound1 = as<int>(Rb1);
            break;
        }
        default: {
            stop("bound1 must be of type numeric or integer");
        }
    }
    if (bound1 < 1 || bound1 > pow(2,31) - 1) {stop("bound1 must be a positive number less than 2^31");}

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
                bound2 = as<int>(Rb2);
                break;
            }
            case INTSXP: {
                bound2 = as<int>(Rb2);
                break;
            }
            default: {
                stop("bound2 must be of type numeric or integer");
            }
        }
        if (bound2 < 1 || bound2 > pow(2,31) - 1) {stop("bound2 must be a positive number less than 2^31");}

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
