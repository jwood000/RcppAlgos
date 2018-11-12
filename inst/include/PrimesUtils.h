#ifndef PRIMES_UTILS_H
#define PRIMES_UTILS_H

#include <PrimesSegSieve.h>
#include <cmath>
#include <vector>
#include <algorithm>

// This is the largest multiple of 2*3*5*7 = 210
// that is less than 2^15 = 32768 = 32KB. This
// is the typical size of most CPU's L1 cache
const int L1CacheSize = 32760;
const unsigned long int wheelSize = 48;
const int_fast64_t sz210 = 210;
constexpr unsigned long int sz420 = 2 * sz210;
constexpr int_fast64_t segmentSize = (int_fast64_t) L1CacheSize;
constexpr unsigned long int nWheelsPerSeg = (L1CacheSize / sz210);

static const int_fast64_t wheel210[wheelSize] = {
    10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4,
    2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2};

static const char check210[sz210] = {
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
    0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
    0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
    1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
    0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
    0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
    0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1};

static const unsigned long int remainder210[sz420] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94,
    95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113,
    114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131,
    132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
    150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
    168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185,
    186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203,
    204, 205, 206, 207, 208, 209, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86,
    87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107,
    108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
    126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161,
    162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
    180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
    198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209};

const unsigned long int MaxBucketSize = 1024;

struct Sieve2dData {
    int_fast64_t svPri;
    int_fast64_t nextStrt;
};

struct Bucket {
    Sieve2dData sieve2dPrimes[MaxBucketSize];
    std::size_t bucketSize = 0;
};

constexpr unsigned long int smlPriBsSize = sizeof(smallPrimeBase) / sizeof(smallPrimeBase[0]);

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

const std::vector<double> percInc = {0.25, 0.116, 0.103, 0.0850, 0.0712,
                                     0.0614, 0.0538, 0.0480, 0.0431, 
                                     0.0392, 0.0360, 0.0332, 0.0309};

const std::vector<double> cutPoints = {40000.0, 120000.0, 1000000.0, 10000000.0,
                                       100000000.0, 1000000000.0, 10000000000.0,
                                       100000000000.0, 1000000000000.0,
                                       10000000000000.0, 100000000000000.0,
                                       1000000000000000.0, 10000000000000000.0};

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
