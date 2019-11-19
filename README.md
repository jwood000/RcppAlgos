
# RcppAlgos <img src='/inst/figures/RcppAlgos-logo.png' width="181px" align="right" />

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
  * E.g. Finding all combinations of a vector such that the mean is between two values.
* A highly efficient generalized partition function is employed when `constraintFun = "sum"` and `comparisonFun = "=="` (see examples below).
* Easily generate random samples of combinations/permutations with **combo/permuteSample**.
* Produce results in parallel using the `Parallel`  or `nThreads` arguments. You can also apply each of the five compiled functions given by the argument `constraintFun` in parallel as well.
* GMP support allows for exploration of combinations/permutations of vectors with many elements.
* Offers a variety of high performance number theoretic functions that are useful for problems common in computational mathematics (E.g. **primeSieve**, **primeFactorize**, and **divisorsRcpp** to name a few).

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
system.time(primeSieve(1e15, 1e15 + 1e9, nThreads = 8))
   user  system elapsed 
  4.872   0.807   0.917

  
## Count primes quickly
system.time(print(primeCount(1e14, nThreads = 8)))
[1] 3204941750802
   user  system elapsed 
 50.121   0.054   6.688


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


## Get combinations such that the product is between
## 3600 and 4000 (including 3600 but not 4000)
comboGeneral(5, 7, TRUE, constraintFun = "prod",
             comparisonFun = c(">=","<"),
             limitConstraints = c(3600, 4000),
             keepResults = TRUE)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
[1,]    1    2    3    5    5    5    5 3750
[2,]    1    3    3    4    4    5    5 3600
[3,]    1    3    4    4    4    4    5 3840
[4,]    2    2    3    3    4    5    5 3600
[5,]    2    2    3    4    4    4    5 3840
[6,]    3    3    3    3    3    3    5 3645
[7,]    3    3    3    3    3    4    4 3888


## Find all combinations of the vector c(121, 126, ..., 216, 221)
## of length 13 such that the sum is equal 2613
system.time(parts <- comboGeneral(seq(121L, 221L, 5L), 13, TRUE,
                                  constraintFun = "sum", 
                                  comparisonFun = "==", 
                                  limitConstraints = 2613))
   user  system elapsed 
  0.008   0.001   0.009
  
head(parts)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
[1,]  121  121  161  221  221  221  221  221  221   221   221   221   221
[2,]  121  121  166  216  221  221  221  221  221   221   221   221   221
[3,]  121  121  171  211  221  221  221  221  221   221   221   221   221
[4,]  121  121  171  216  216  221  221  221  221   221   221   221   221
[5,]  121  121  176  206  221  221  221  221  221   221   221   221   221
[6,]  121  121  176  211  216  221  221  221  221   221   221   221   221

dim(parts)
[1] 119546     13

## Over 500 million possible results
prettyNum(comboCount(seq(121L, 221L, 5L), 13, TRUE), big.mark = ",")
[1] "573,166,440"
```

## Getting Started

* Combinatorial Sampling
* Computational Mathematics Overview
* Combinations and Permutations under Constraints
* Combination and Permutation Basics

## Contact

I welcome any and all feedback. If you would like to report a bug, have a question, or have suggestions for possible improvements, please file an [issue](<https://github.com/jwood000/RcppAlgos/issues>).
