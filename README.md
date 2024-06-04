
# RcppAlgos <img src='man/figures/RcppAlgos-logo.png' width="150px" align="right" />

<!-- badges: start -->
[![CRAN status](<https://www.r-pkg.org/badges/version/RcppAlgos>)](<https://cran.r-project.org/package=RcppAlgos>)
![](<http://cranlogs.r-pkg.org/badges/RcppAlgos?color=orange>)
![](<http://cranlogs.r-pkg.org/badges/grand-total/RcppAlgos?color=brightgreen>)
[![Codacy Badge](<https://app.codacy.com/project/badge/Grade/e7fef773f6514aa4a2decda9adf57ae8>)](<https://app.codacy.com/gh/jwood000/RcppAlgos/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade>)
[![Dependencies](<https://tinyverse.netlify.app/badge/RcppAlgos>)](<https://cran.r-project.org/package=RcppAlgos>)
[![R-CMD-check](https://github.com/jwood000/RcppAlgos/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/jwood000/RcppAlgos/actions/workflows/R-CMD-check.yml)
[![Codecov test coverage](https://codecov.io/gh/jwood000/RcppAlgos/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jwood000/RcppAlgos?branch=main)
<!-- badges: end -->

A collection of high performance functions and iterators implemented in C++ for solving problems in combinatorics and computational mathematics.

## Featured Functions

  - **`{combo|permute}General`**: Generate all combinations/permutations of a vector (including [multisets](<https://en.wikipedia.org/wiki/Multiset>)) meeting specific criteria.
  - **`{partitions|compositions}General`**: Efficient algorithms for partitioning numbers under various constraints
  - **`{combo|permute|partitions|compositions}Sample`**: Generate reproducible random samples
  - **`{combo|permute|partitions|compositions}Iter`**: Flexible iterators allow for bidirectional iteration as well as random access.
  - **`primeSieve`**: Fast prime number generator
  - **`primeCount`**: Prime counting function using [Legendre's formula](<http://mathworld.wolfram.com/LegendresFormula.html>)

The `primeSieve` function and the `primeCount` function are both based off of the excellent work by [Kim Walisch](<https://github.com/kimwalisch>). The respective repos can be found here: [kimwalisch/primesieve](<https://github.com/kimwalisch/primesieve>); [kimwalisch/primecount](<https://github.com/kimwalisch/primecount>)

Additionally, many of the sieving functions make use of the fast integer division library [libdivide](<https://github.com/ridiculousfish/libdivide>) by [ridiculousfish](<https://github.com/ridiculousfish>).

## Benchmarks

* [High Performance Benchmarks](<https://jwood000.github.io/RcppAlgos/articles/HighPerformanceBenchmarks.html>)

## Installation

``` r
install.packages("RcppAlgos")

## install the development version
devtools::install_github("jwood000/RcppAlgos")
```

## Basic Usage

### Combinatorics

``` r
## Find all 3-tuples combinations of 1:4
comboGeneral(4, 3)
#>      [,1] [,2] [,3]
#> [1,]   1    2    3
#> [2,]   1    2    4
#> [3,]   1    3    4
#> [4,]   2    3    4


## Alternatively, iterate over combinations
a = comboIter(4, 3)
a@nextIter()
#> [1] 1 2 3

a@back()
#> [1] 2 3 4

a[[2]]
#> [1] 1 2 4


## Pass any atomic type vector
permuteGeneral(letters, 3, upper = 4)
#>      [,1] [,2] [,3]
#> [1,] "a"  "b"  "c"
#> [2,] "a"  "b"  "d"
#> [3,] "a"  "b"  "e"
#> [4,] "a"  "b"  "f"


## Flexible partitioning algorithms
partitionsGeneral(0:5, 3, freqs = rep(1:2, 3), target = 6)
#>      [,1] [,2] [,3]
#> [1,]    0    1    5
#> [2,]    0    2    4
#> [3,]    0    3    3
#> [4,]    1    1    4
#> [5,]    1    2    3


## And compositions
compositionsGeneral(0:3, repetition = TRUE)
#>      [,1] [,2] [,3]
#> [1,]    0    0    3
#> [2,]    0    1    2
#> [3,]    0    2    1
#> [4,]    1    1    1


## Generate a reproducible sample
comboSample(10, 8, TRUE, n = 5, seed = 84)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    3    3    3    6    6   10   10   10
#> [2,]    1    3    3    4    4    7    9   10
#> [3,]    3    7    7    7    9   10   10   10
#> [4,]    3    3    3    9   10   10   10   10
#> [5,]    1    2    2    3    3    4    4    7


## Get combinations such that the product is between
## 3600 and 4000 (including 3600 but not 4000)
comboGeneral(5, 7, TRUE, constraintFun = "prod",
             comparisonFun = c(">=","<"),
             limitConstraints = c(3600, 4000),
             keepResults = TRUE)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    1    2    3    5    5    5    5 3750
#> [2,]    1    3    3    4    4    5    5 3600
#> [3,]    1    3    4    4    4    4    5 3840
#> [4,]    2    2    3    3    4    5    5 3600
#> [5,]    2    2    3    4    4    4    5 3840
#> [6,]    3    3    3    3    3    3    5 3645
#> [7,]    3    3    3    3    3    4    4 3888


## We can even iterate over constrained cases. These are
## great when we don't know how many results there are upfront.
## Save on memory and still at the speed of C++!!
p = permuteIter(5, 7, TRUE, constraintFun = "prod",
                comparisonFun = c(">=","<"),
                limitConstraints = c(3600, 4000),
                keepResults = TRUE)

## Get the next n results
t <- p@nextNIter(1048)

## N.B. keepResults = TRUE adds the 8th column
tail(t)
#>         [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1043,]    5    4    4    3    4    1    4 3840
#> [1044,]    5    4    4    3    4    4    1 3840
#> [1045,]    5    4    4    4    1    3    4 3840
#> [1046,]    5    4    4    4    1    4    3 3840
#> [1047,]    5    4    4    4    3    1    4 3840
#> [1048,]    5    4    4    4    3    4    1 3840

## Continue iterating from where we left off
p@nextIter()
#> [1]    5    4    4    4    4    1    3 3840

p@nextIter()
#> [1]    5    4    4    4    4    3    1 3840

p@nextIter()
#> [1]    2    2    3    3    4    5    5 3600

## N.B. totalResults and totalRemaining are NA because there is no
## closed form solution for determining this.
p@summary()
#> $description
#> [1] "Permutations with repetition of 5 choose 7 where the prod is between 3600 and 4000"
#> 
#> $currentIndex
#> [1] 1051
#> 
#> $totalResults
#> [1] NA
#> 
#> $totalRemaining
#> [1] NA
```

### Computational Mathematics

```r
## Generate prime numbers
primeSieve(50)
#> [1]  2  3  5  7 11 13 17 19 23 29 31 37 41 43 47

## Many of the functions can produce results in
## parallel for even greater performance
p = primeSieve(1e15, 1e15 + 1e8, nThreads = 4)

head(p)
#> [1] 1000000000000037 1000000000000091 1000000000000159
#> [4] 1000000000000187 1000000000000223 1000000000000241
tail(p)
#> [1] 1000000099999847 1000000099999867 1000000099999907
#> [4] 1000000099999919 1000000099999931 1000000099999963


## Count prime numbers less than n
primeCount(1e10)
#> [1] 455052511

## Get the prime factorization
set.seed(24028)
primeFactorize(sample(1e15, 3), namedList = TRUE)
#> $`701030825091514`
#> [1]             2           149 2352452433193
#> 
#> $`83054168594779`
#> [1]  3098071 26808349
#> 
#> $`397803024735610`
#> [1]            2            5           13           13 235386405169
```

## Further Reading

* [Function Documentation](<https://jwood000.github.io/RcppAlgos/reference/index.html>)
* [Computational Mathematics Overview](<https://jwood000.github.io/RcppAlgos/articles/ComputationalMathematics.html>)
* [Combination and Permutation Basics](<https://jwood000.github.io/RcppAlgos/articles/GeneralCombinatorics.html>)
* [Combinatorial Sampling](<https://jwood000.github.io/RcppAlgos/articles/CombinatorialSampling.html>)
* [Constraints, Partitions, and Compositions](<https://jwood000.github.io/RcppAlgos/articles/CombPermConstraints.html>)
* [Attacking Problems Related to the Subset Sum Problem](<https://jwood000.github.io/RcppAlgos/articles/SubsetSum.html>)
* [Combinatorial Iterators in RcppAlgos](<https://jwood000.github.io/RcppAlgos/articles/CombinatoricsIterators.html>)
* [Cartesian Products and Partitions of Groups](<https://jwood000.github.io/RcppAlgos/articles/OtherCombinatorics.html>)

## Why **`RcppAlgos`** but no **`Rcpp`**?

Previous versions of `RcppAlgos` relied on [Rcpp](<https://github.com/RcppCore/Rcpp>) to ease the burden of exposing C++ to R. While the current version of `RcppAlgos` does not utilize `Rcpp`, it would not be possible without the myriad of excellent contributions to `Rcpp`.

## Contact

If you would like to report a bug, have a question, or have suggestions for possible improvements, please file an [issue](<https://github.com/jwood000/RcppAlgos/issues>).
