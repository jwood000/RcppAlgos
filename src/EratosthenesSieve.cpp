#include "EratosthenesSieve.h"
#include "CleanConvert.h"

// [[Rcpp::export]]
SEXP EratosthenesRcpp(SEXP Rb1, SEXP Rb2, SEXP RNumThreads, int maxCores, int maxThreads) {
    
    double bound1, bound2;
    int_fast64_t myMax, myMin;
    CleanConvert::convertPrimitive(Rb1, bound1, "bound1 must be of type numeric or integer", false);
    
    if (bound1 <= 0 || bound1 > PrimeSieve::Significand53)
        Rcpp::stop("bound1 must be a positive number less than 2^53");
    
    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, "bound2 must be of type numeric or integer", false);
    }
    
    if (bound2 <= 0 || bound2 > PrimeSieve::Significand53)
        Rcpp::stop("bound2 must be a positive number less than 2^53");
    
    if (bound1 > bound2) {
        myMax = static_cast<int_fast64_t>(bound1);
        myMin = static_cast<int_fast64_t>(bound2);
    } else {
        myMax = static_cast<int_fast64_t>(bound2);
        myMin = static_cast<int_fast64_t>(bound1);
    }
    
    myMin = std::ceil(myMin);
    myMax = std::floor(myMax);
        
    if (myMax <= 1)
        return Rcpp::IntegerVector();
    
    if (myMin <= 2) myMin = 1;
    if (myMin == myMax) {++myMax;}
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads must be of type numeric or integer");
    
    std::size_t numPrimes = 0u;
    std::vector<unsigned long int> runningCount;
    runningCount.push_back(0u);
    unsigned long int numSects = nThreads;
    bool Parallel = false;
    
    if (myMax > INT_MAX) {
        std::vector<std::vector<double>> primeList(numSects, std::vector<double>());
        std::vector<double> tempPrime;
        
        PrimeSieve::PrimeMaster(myMin, myMax, tempPrime, primeList, 
                                Parallel, nThreads, maxThreads, maxCores);
        
        if (Parallel) {
            for (int i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }
            
            Rcpp::NumericVector primes(numPrimes);
            Rcpp::NumericVector::iterator priBeg = primes.begin();
            
            for (int i = 0; i < numSects; ++i)
                std::copy(primeList[i].begin(), primeList[i].end(), priBeg + runningCount[i]);
            
            return primes;
        } else {
            return Rcpp::wrap(tempPrime);
        }
    } else {
        std::vector<std::vector<int_fast32_t>> primeList(numSects, std::vector<int_fast32_t>());
        std::vector<int_fast32_t> tempPrime;
        
        PrimeSieve::PrimeMaster(myMin, myMax, tempPrime, primeList, 
                                Parallel, nThreads, maxThreads, maxCores);
        
        if (Parallel) {
            for (int i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }
            
            Rcpp::IntegerVector primes(numPrimes);
            Rcpp::IntegerVector::iterator priBeg = primes.begin();
            
            for (int i = 0; i < numSects; ++i)
                std::copy(primeList[i].begin(), primeList[i].end(), priBeg + runningCount[i]);
            
            return primes;
        } else {
            return Rcpp::wrap(tempPrime);
        }
    }
}
