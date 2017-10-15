# algosRcpp

Overview
---------
A collection of optimized functions implemented in C++ with Rcpp for solving problems in combinatorics and computational mathematics.

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
comboGeneral(5, 3, v = c(2, 3, 5, 7, 11), 
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
```
