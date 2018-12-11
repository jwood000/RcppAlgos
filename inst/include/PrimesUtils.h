#ifndef PRIMES_UTILS_H
#define PRIMES_UTILS_H

#include <PrimesSegSieve.h>
#include <Rcpp.h>
#include <cmath>

// This is the largest multiple of 2*3*5*7 = 210
// that is less than 2^15 = 32768 = 32KB. This
// is the typical size of most CPU's L1 cache
constexpr int L1CacheSize = 32760;
constexpr unsigned long int wheelSize = 48;
constexpr std::size_t nWheelsPerSeg = static_cast<std::size_t>(L1CacheSize / 210);

static const int_fast64_t wheel210[wheelSize] = {
    10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4,
    2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2};

int_fast64_t smallCut = 1000000000;

const int_fast64_t maxPriPer = 200000000;    // 2e8
// maxPriPer was obtained by noting that the number pi(2e8) ~= 11e6.
// Empirical tests with std::vector and push_back shows significant
// performance decrease with larger vectors. E.g. 
//
//    void push_backTest(NumericVector x, bool verbose = false) {
//       std::vector<double> v;
//        v.reserve(x.size());
//        
//        for (std::size_t i = 0; i < x.size(); ++i)
//            v.push_back(x[i]);
//    }
//
//                            ## millions
//    v1e6 = sample(1e8, 1e6);                 v2e6 = sample(1e8, 2e6)  
//    system.time(push_backTest(v2e6)) ~= 2 * system.time(push_backTest(v1e6))
//
// And the same is true for smaller cases (i.e. we get linear performance), however
// when we get to around 15e6, the performance drops significantly
//
//                            ## ten-millions
//    v1e7 = sample(1e8, 1e7);                  v2e7 = sample(1e8, 2e7)  
//    system.time(push_backTest(v2)) ~= 4 * system.time(push_backTest(v1))
//
// This seems to be due to memory performance. We see the same decrease in
// performance when we use IntegerVector (coupled with std::vector<int>) 
// instead of NumericVector (and std::vector<double>) only at a point that
// is around twice as large... i.e. The performance is linear until about
// 30e6 (it was 15e6 before). Note that a double takes up 8 bytes whereas
// and integer only takes 4 bytes. With this in mind, the const we have
// chosen may not be a good fit on all machines. This topic needs revisting.

// These numbers were obtained empirically using the prime number theorem
// along with a prime counting function. Values were computed in each
// range below in cutPoints. For example in the range [40000, 120000),
// if we increment by 10 and calculate the following: 
//                       pi(x) - x / log(x) = e
// We then calculate the percentage of x / log(x) or : e / (x / log(x))
// for every value and find the max:
//             a <- sapply(seq(40000, 120000, 10), primeCount)
//             b <- sapply(seq(40000, 120000, 10), function(x) x/log(x))
//             myDiff <- a - b;  max(myDiff / b) ## [1] 0.1153694

const std::vector<double> percInc = {0.2500, 0.1160, 0.1030, 0.0850, 0.0712,
                                     0.0614, 0.0538, 0.0495, 0.0480, 0.0431, 
                                     0.0392, 0.0360, 0.0332, 0.0309, 0.0292};

const std::vector<double> cutPoints = {            40000.0,           120000.0,
                                                 1000000.0,         10000000.0,
                                               100000000.0,       1000000000.0,
                                              5000000000.0,      10000000000.0,
                                            100000000000.0,    1000000000000.0,
                                          10000000000000.0,  100000000000000.0,
                                        1000000000000000.0, 5000000000000000.0,
                                       10000000000000000.0};

// The following function is based off of the prime number theorem
std::size_t EstimatePiPrime(double minNum, double maxNum) {
    std::vector<double>::const_iterator it = std::upper_bound(cutPoints.begin(),
                                                              cutPoints.end(),
                                                              maxNum);
    std::size_t myIndex = it - cutPoints.begin();
    double dblRes = std::ceil((maxNum / log(maxNum)) * (1 + percInc[myIndex]));
    
    if (minNum > 1000)
        dblRes -= std::floor((minNum / log(minNum)) * (1 + percInc[myIndex]));
    
    std::size_t result = dblRes;
    return result;
}

#endif
