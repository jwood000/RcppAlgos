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
## Find all 3-tuples combinations without 
## repetition of the numbers c(1, 2, 3, 4).
comboGeneral(4, 3)
     [,1] [,2] [,3]
[1,]   1    2    3
[2,]   1    2    4
[3,]   1    3    4
[4,]   2    3    4


## Find all 3-tuples permutations without
## repetition of the numbers c(1, 2, 3, 4).
permuteGeneral(4,3)
      [,1] [,2] [,3]
 [1,]    1    2    3
 [2,]    2    1    3
 [3,]    3    1    2
 [4,]    1    3    2
  .      .    .    .
  .      .    .    .
[21,]    4    2    3
[22,]    2    4    3
[23,]    3    4    2
[24,]    4    3    2

## Find all 3-tuples combinations without repetition of
## c(2, 3, 5, 7, 11), such that the product is less than 130.
## The last column is the resulting sum of each combination
comboGeneral(v = c(2, 3, 5, 7, 11),
             m = 3, 
             constraintFun = "prod", 
             comparisonFun = "<", 
             limitConstraints = 130,
             keepResults = TRUE)
     [,1] [,2] [,3] [,4]
[1,]    2    3    5   30
[2,]    2    3    7   42
[3,]    2    3   11   66
[4,]    2    5    7   70
[5,]    2    5   11  110
[6,]    3    5    7  105


## Generate some random data
set.seed(101)
s <- sample(500, 20)
## Find all 5-tuples permutations without repetition
## of s (defined above) such that the sum is equal to 1176.
p <- permuteGeneral(v = s,
                    m = 5, 
                    constraintFun = "sum", 
                    comparisonFun = "==", 
                    limitConstraints = 1176,
                    keepResults = TRUE)
                    
head(p)
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,]   19   22  327  354  454 1176
[2,]   22   19  327  354  454 1176
[3,]  327   19   22  354  454 1176
[4,]   19  327   22  354  454 1176
[5,]   22  327   19  354  454 1176
[6,]  327   22   19  354  454 1176
  
tail(p)
        [,1] [,2] [,3] [,4] [,5] [,6]
[3955,]  199  187  354  287  149 1176
[3956,]  187  199  354  287  149 1176
[3957,]  354  199  187  287  149 1176
[3958,]  199  354  187  287  149 1176
[3959,]  187  354  199  287  149 1176
[3960,]  354  187  199  287  149 1176


## Generate all permutations of a vector with specific
## length of repetition for each element
permuteGeneral(3, freqs = c(1,2,2))
      [,1] [,2] [,3] [,4] [,5]
 [1,]    1    2    2    3    3
 [2,]    1    2    3    2    3
 [3,]    1    2    3    3    2
 [4,]    1    3    2    2    3
 [5,]    1    3    2    3    2
  .      .    .    .    .    .
  .      .    .    .    .    .
[26,]    3    2    3    1    2
[27,]    3    2    3    2    1
[28,]    3    3    1    2    2
[29,]    3    3    2    1    2
[30,]    3    3    2    2    1


## All combinatoric functions work with factors
facPerms <- permuteGeneral(factor(c("low", "med", "high")),
                                     freqs = c(1,4,2))
facPerms[1:5, ]
     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,] high low  low  low  low  med  med 
[2,] high low  low  low  med  low  med 
[3,] high low  low  low  med  med  low 
[4,] high low  low  med  low  low  med 
[5,] high low  low  med  low  med  low 
Levels: high low med


## They are very efficient as well!!
system.time(permuteGeneral(4, freqs = 5:2))
## That's over 2.5 million permutations instantly!!


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
