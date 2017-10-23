# RcppAlgos

Overview
---------
A collection of optimized functions implemented in C++ with Rcpp for solving problems in combinatorics and computational mathematics.

Installation
------------

``` r
# install.packages("devtools")
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
system.time(myPs <- primeSieve(1000001000, 10^9))

myPs
 [1] 1000000007 1000000009 1000000021 1000000033 1000000087
 [6] 1000000093 1000000097 1000000103 1000000123 1000000181
[11] 1000000207 1000000223 1000000241 1000000271 1000000289
[16] 1000000297 1000000321 1000000349 1000000363 1000000403
[21] 1000000409 1000000411 1000000427 1000000433 1000000439
[26] 1000000447 1000000453 1000000459 1000000483 1000000513
[31] 1000000531 1000000579 1000000607 1000000613 1000000637
[36] 1000000663 1000000711 1000000753 1000000787 1000000801
[41] 1000000829 1000000861 1000000871 1000000891 1000000901
[46] 1000000919 1000000931 1000000933 1000000993

## Object created is small
object.size(myPs)
240 bytes

## Generate number of coprime elements for many numbers
myPhis <- eulerPhiSieve(20)
names(myPhis) <- 1:20
myPhis
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
 1  1  2  2  4  2  6  4  6  4 10  4 12  6  8  8 16  6 18  8
```
