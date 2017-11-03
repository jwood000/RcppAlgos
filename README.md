# RcppAlgos

Overview
---------
A collection of optimized functions implemented in C++ with Rcpp for solving problems in combinatorics and computational mathematics.

Installation
------------

``` r
install.packages("RcppAlgos")

## Or install the development version
devtools::install_github("jwood000/RcppAlgos")
```

Usage
-----
``` r
## Find all 3-tuples without repetition of the numbers c(1, 2, 3, 4).
comboGeneral(4, 3)
     [,1] [,2] [,3]
[1,]   1    2    3
[2,]   1    2    4
[3,]   1    3    4
[4,]   2    3    4


## Find all 3-tuples without repetition of c(2, 3, 5, 7, 11), such that the product is less than 130.
comboGeneral(v = c(2, 3, 5, 7, 11),
             m = 3, 
             constraintFun = "prod", 
             comparisonFun = "<", 
             limitConstraints = 130)
     [,1] [,2] [,3]
[1,]    2    3    5
[2,]    2    3    7
[3,]    2    3   11
[4,]    2    5    7
[5,]    2    5   11
[6,]    3    5    7


 ## get prime factorization for every number from 1 to n
 primeFactorizationList(5)
[[1]]
integer(0)

[[2]]
[1] 2

[[3]]
[1] 3

[[4]]
[1] 2 2

[[5]]
[1] 5

## Quickly generate large primes over small interval
options(scipen = 50)
system.time(myPs <- primeSieve(10^13+10^3, 10^13))

myPs
 [1] 10000000000037 10000000000051 10000000000099 10000000000129
 [5] 10000000000183 10000000000259 10000000000267 10000000000273
 [9] 10000000000279 10000000000283 10000000000313 10000000000343
[13] 10000000000391 10000000000411 10000000000433 10000000000453
[17] 10000000000591 10000000000609 10000000000643 10000000000649
[21] 10000000000657 10000000000687 10000000000691 10000000000717
[25] 10000000000729 10000000000751 10000000000759 10000000000777
[29] 10000000000853 10000000000883 10000000000943 10000000000957
[33] 10000000000987 10000000000993

## Object created is small
object.size(myPs)
312 bytes

## Generate number of coprime elements for many numbers
myPhis <- eulerPhiSieve(20)
names(myPhis) <- 1:20
myPhis
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
 1  1  2  2  4  2  6  4  6  4 10  4 12  6  8  8 16  6 18  8
```
