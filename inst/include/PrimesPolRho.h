#ifndef PRIMES_POL_RHO_H
#define PRIMES_POL_RHO_H

#include <array>

// The first 550 prime numbers (i.e. the first prime numbers less than 4000).
// Starting with the first prime number 2, add each successive element in
// primesDiffPR to obtain the next prime number (e.g. p = 2 + primesDiffPR[0] = 3,
// p = p + primesDiffPR[1] = 5, etc.)
// 
// N.B. The last prime in primesDiffPR is 3989, Thus the first prime omitted
// is 4001. That is:
//
// 3989 = std::accumulate(primesDiffPR.cbegin(), primesDiffPR.cend(), 2);
//                        - and -
// gmp::nextprime(3989) -->> 4001

constexpr int64_t FirstOmittedPrime = 4001;

static const std::array<int, 549> primesDiffPR = {{
    1,2,2,4,2,4,2,4,6,2,6,4,2,4,6,6,2,6,4,2,6,4,6,8,4,2,4,
    2,4,14,4,6,2,10,2,6,6,4,6,6,2,10,2,4,2,12,12,4,2,4,6,
    2,10,6,6,6,2,6,4,2,10,14,4,2,4,14,6,10,2,4,6,8,6,6,4,
    6,8,4,8,10,2,10,2,6,4,6,8,4,2,4,12,8,4,8,4,6,12,2,18,
    6,10,6,6,2,6,10,6,6,2,6,6,4,2,12,10,2,4,6,6,2,12,4,6,
    8,10,8,10,8,6,6,4,8,6,4,8,4,14,10,12,2,10,2,4,2,10,14,
    4,2,4,14,4,2,4,20,4,8,10,8,4,6,6,14,4,6,6,8,6,12,4,6,
    2,10,2,6,10,2,10,2,6,18,4,2,4,6,6,8,6,6,22,2,10,8,10,
    6,6,8,12,4,6,6,2,6,12,10,18,2,4,6,2,6,4,2,4,12,2,6,34,
    6,6,8,18,10,14,4,2,4,6,8,4,2,6,12,10,2,4,2,4,6,12,12,
    8,12,6,4,6,8,4,8,4,14,4,6,2,4,6,2,6,10,20,6,4,2,24,4,
    2,10,12,2,10,8,6,6,6,18,6,4,2,12,10,12,8,16,14,6,4,2,
    4,2,10,12,6,6,18,2,16,2,22,6,8,6,4,2,4,8,6,10,2,10,14,
    10,6,12,2,4,2,10,12,2,16,2,6,4,2,10,8,18,24,4,6,8,16,
    2,4,8,16,2,4,8,6,6,4,12,2,22,6,2,6,4,6,14,6,4,2,6,4,6,
    12,6,6,14,4,6,12,8,6,4,26,18,10,8,4,6,2,6,22,12,2,16,
    8,4,12,14,10,2,4,8,6,6,4,2,4,6,8,4,2,6,10,2,10,8,4,14,
    10,12,2,6,4,2,16,14,4,6,8,6,4,18,8,10,6,6,8,10,12,14,4,
    6,6,2,28,2,10,8,4,14,4,8,12,6,12,4,6,20,10,2,16,26,4,2,
    12,6,4,12,6,8,4,8,22,2,4,2,12,28,2,6,6,6,4,6,2,12,4,12,
    2,10,2,16,2,16,6,20,16,8,4,2,4,2,22,8,12,6,10,2,4,6,2,
    6,10,2,12,10,2,10,14,6,4,6,8,6,6,16,12,2,4,14,6,4,8,10,
    8,6,6,22,6,2,10,14,4,6,18,2,10,14,4,2,10,14,4,8,18,4,6,
    2,4,6,2,12,4,20,22}};

#endif
