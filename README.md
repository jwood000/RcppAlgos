
# RcppAlgos <img src="./inst/figures/RcppAlgos-logo.png" width="181px" align="right" />

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/jwood000/RcppAlgos.svg?branch=master)](https://travis-ci.com/jwood000/RcppAlgos)
![](http://cranlogs.r-pkg.org/badges/RcppAlgos?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/RcppAlgos?color=brightgreen)
[![Coverage status](https://codecov.io/gh/jwood000/RcppAlgos/branch/master/graph/badge.svg)](https://codecov.io/github/jwood000/RcppAlgos?branch=master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/47c85553b9bd4479a722857c1f2fdcf3)](https://app.codacy.com/app/jwood000/RcppAlgos?utm_source=github.com&utm_medium=referral&utm_content=jwood000/RcppAlgos&utm_campaign=Badge_Grade_Settings)
[![Dependencies](https://tinyverse.netlify.com/badge/RcppAlgos)](https://cran.r-project.org/package=RcppAlgos)
<!-- badges: end -->

A collection of high performance functions implemented in C++ with Rcpp for solving problems in combinatorics and computational mathematics. Utilizes the [RcppThread](<https://github.com/tnagler/RcppThread>) library for easy access to thread safe multithreading. We also make use of the [RMatrix.h](<https://github.com/RcppCore/RcppParallel/blob/master/inst/include/RcppParallel/RMatrix.h>) header file from [RcppParallel](<https://github.com/RcppCore/RcppParallel>) for thread safe accessors for Rcpp matrices.

## Features

* Generate all combinations/permutations of a vector (including [multisets](<https://en.wikipedia.org/wiki/Multiset>)) meeting specific criteria with **combo/permuteGeneral**.
* A highly efficient generalized partition function is employed when `constraintFun = "sum"` and `comparisonFun = "=="`.
* Produce results in parallel using the `Parallel` argument. You can also apply each of the five compiled functions given by the argument `constraintFun` in parallel as well. E.g. Obtaining the row sums of all combinations
* GMP support allows for exploration of combinations/permutations of vectors with many elements.
* Easily generate random samples of combinations/permutations in parallel with **combo/permuteSample**
* Offers a variety high performance number theoretic functions that are useful for problems common in computational mathematics
* Generates all primes less than a **billion in under 0.4 seconds** with **primeSieve**
* Count the primes less than **1e15 in under 50 seconds** with **primeCount**

The `primeSieve` function and the `primeCount` function are both based off of the excellent work by [Kim Walisch](<https://github.com/kimwalisch>). The respective repos can be found here: [kimwalisch/primesieve](<https://github.com/kimwalisch/primesieve>); [kimwalisch/primecount](<https://github.com/kimwalisch/primecount>)

Additionally, many of the sieving functions make use of the fast integer division library [libdivide](<https://github.com/ridiculousfish/libdivide>) by [ridiculousfish](<https://github.com/ridiculousfish>).

## Installation

``` r
install.packages("RcppAlgos")

## Or install the development version
devtools::install_github("jwood000/RcppAlgos")
```

## Usage

``` r
## Generate primes in a range quickly
system.time(b <- primeSieve(1e15, 1e15 + 1e9, nThreads = 8))
   user  system elapsed 
  4.898   0.889   0.953

  
## Count primes quickly
system.time(print(primeCount(1e14, nThreads = 8)))
[1] 3204941750802
   user  system elapsed 
 53.371   0.117   7.221


## Completely factor a vector of numbers 
set.seed(42)
myNums <- sample(1e12, 3)
allDivs <- divisorsRcpp(myNums, namedList = TRUE)


## Find all 3-tuples combinations without 
## repetition of the numbers c(1, 2, 3, 4).
comboGeneral(4, 3)
     [,1] [,2] [,3]
[1,]   1    2    3
[2,]   1    2    4
[3,]   1    3    4
[4,]   2    3    4


## Find the first four 3-tuples permutations without
## repetition of the numbers c(1, 2, 3, 4).
permuteGeneral(4, 3, upper = 4)
      [,1] [,2] [,3]
 [1,]    1    2    3
 [2,]    1    2    4
 [3,]    1    3    2
 [4,]    1    3    4


## Use lower to specify which lexicographical combination
## to start generating from (if upper not specified, we 
## go to the last result i.e. comboCount(4, 3, T) = 20)
comboGeneral(4, 3, TRUE, lower = 17)
      [,1] [,2] [,3]
[1,]    3    3    3
[2,]    3    3    4
[3,]    3    4    4
[4,]    4    4    4
  

## They are very efficient
system.time(comboGeneral(25, 13))
   user  system elapsed 
  0.104   0.054   0.158

## Even faster in parallel
system.time(comboGeneral(25, 13, nThreads = 8))
   user  system elapsed 
  0.166   0.220   0.055


## Generate a reproducible sample of combinations
comboSample(10, 8, TRUE, n = 5, seed = 84)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
[1,]    3    3    3    6    6   10   10   10
[2,]    1    3    3    4    4    7    9   10
[3,]    3    7    7    7    9   10   10   10
[4,]    3    3    3    9   10   10   10   10
[5,]    1    2    2    3    3    4    4    7
```

## Getting Started

* [Combinatorial Sampling]
* [Computational Mathematics Overview]
* [Combinations and Permutations under Constraints]
* [Combination and Permutation Basics]

## Contact

I welcome any and all feedback. If you would like to report a bug, have a question, or have suggestions for possible improvements, please contact me here: jwood000@gmail.com
