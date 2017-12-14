[![Rdoc](http://www.rdocumentation.org/badges/version/RcppAlgos)](http://www.rdocumentation.org/packages/RcppAlgos)
![](http://cranlogs.r-pkg.org/badges/RcppAlgos?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/RcppAlgos?color=brightgreen)

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

### Common Combinatorial Functions
Easily executed with a very simple interface.
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
permuteGeneral(4, 3)
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

## For combinations/permutations with repetition, simply
## set the repetition argument to TRUE
comboGeneral(4, 3, repetition = TRUE)
      [,1] [,2] [,3]
 [1,]    1    1    1
 [2,]    1    1    2
 [3,]    1    1    3
  .      .    .    .
  .      .    .    .
  
## They are very efficient
system.time(comboGeneral(25,13))
   user  system elapsed 
  0.174   0.060   0.232 

nrow(comboGeneral(25,13))
[1] 5200300
```

### Using Constraints
Oftentimes, one needs to find combinations/permutations that meet certain requirements.

``` r
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


## Maybe you are curious to see the results of applying a function
## without any constraints. Simply pick the function you wish
## to be applied, and set keepResults to TRUE.
set.seed(99)
mySamp <- rnorm(5, 100, 5)
mySamp
[1] 101.06981 102.39829 100.43914 102.21929  98.18581
comboGeneral(mySamp, m = 4,
             repetition = TRUE,
             constraintFun = "sum",
             keepResults = TRUE)
           [,1]      [,2]      [,3]      [,4]     [,5]
 [1,]  98.18581  98.18581  98.18581  98.18581 392.7432
 [2,]  98.18581  98.18581  98.18581 100.43914 394.9966
 [3,]  98.18581  98.18581  98.18581 101.06981 395.6272
 [4,]  98.18581  98.18581  98.18581 102.21929 396.7767
  .        .         .         .         .         .
  .        .         .         .         .         .
[67,] 102.21929 102.21929 102.21929 102.39829 409.0562
[68,] 102.21929 102.21929 102.39829 102.39829 409.2352
[69,] 102.21929 102.39829 102.39829 102.39829 409.4142
[70,] 102.39829 102.39829 102.39829 102.39829 409.5932
```

### Working with Mulitsets
Sometimes, the standard combination/permutation functions don't quite get us to our desired goals. For
example, one may need all permutations of a vector with some of the elements repeated a specific
amount of times (i.e. a multiset). Consider the following vector `a <- c(1,1,1,1,2,2,2,7,7,7,7,7)` and one
would like to find permutations of `a` of length 6. Using traditional methods, we would need to generate all
permutations, then eliminate duplicate values. Even considering that `permuteGeneral` is very efficient,
this appoach is clunky and not as fast as it could be. Observe:

``` r
getPermsWithSpecificRepetition <- function(z, n) {
    b <- permuteGeneral(z, n)
    myDupes <- duplicated(apply(b, 1, paste, collapse=""))
    b[!myDupes, ]
}

system.time(getPermsWithSpecificRepetition(a, 6))
   user  system elapsed 
  4.300   0.028   4.331
```


#### Enter _freqs_
Situtations like this call for the use of the `freqs` argument. Simply, enter the number
of times each unique element is repeated and Voila!

``` r
system.time(test2 <- permuteGeneral(unique(a), 6, freqs = rle(a)$lengths))
   user  system elapsed 
      0       0       0
      
identical(test[do.call(order,as.data.frame(test)),],
           test2[do.call(order,as.data.frame(test2)),])
[1] TRUE
```
 
 Here are some more general examples with multisets:
 
``` r
## Generate all permutations of a vector with specific
## length of repetition for each element (i.e. multiset)
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

## If you only want a certain length
permuteGeneral(3, 2, freqs = c(1,2,2))
     [,1] [,2]
[1,]    1    2
[2,]    2    1
[3,]    1    3
[4,]    3    1
[5,]    2    2
[6,]    2    3
[7,]    3    2
[8,]    3    3

## or combinations...
comboGeneral(3, 2, freqs = c(1,2,2))
     [,1] [,2]
[1,]    1    2
[2,]    1    3
[3,]    2    2
[4,]    2    3
[5,]    3    3
```

##### All Combinatoric Functions Work with Factors

``` r
facPerms <- permuteGeneral(factor(c("low", "med", "high"),
                                   levels = c("low", "med", "high"),
                                   ordered = TRUE),
                                   freqs = c(1,4,2))
> facPerms[1:5, ]
     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,] low  med  med  med  med  high high
[2,] low  med  med  med  high med  high
[3,] low  med  med  med  high high med 
[4,] low  med  med  high med  med  high
[5,] low  med  med  high med  high med 
Levels: low < med < high
```

Mathematical Computation
-----
`RcppAlgos` comes equipped with several functions for quickly generating essential components
for problems common in computational mathematics.

``` r
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

## If you want the complete factorization from 1 to n, use divisorsList
system.time(allFacs <- divisorsList(10^5))
   user  system elapsed 
  0.059   0.009   0.067
  
names(allFacs) <- 1:10^5

allFacs[c(4339, 15613, 22080, 28372, 40586, 51390, 98248)]
$`4339`
[1]    1 4339

$`15613`
[1]     1    13  1201 15613

$`22080`
 [1]     1     2     3     4     5     6     8    10    12    15
[11]    16    20    23    24    30    32    40    46    48    60
[21]    64    69    80    92    96   115   120   138   160   184
[31]   192   230   240   276   320   345   368   460   480   552
[41]   690   736   920   960  1104  1380  1472  1840  2208  2760
[51]  3680  4416  5520  7360 11040 22080

$`28372`
 [1]     1     2     4    41    82   164   173   346   692  7093
[11] 14186 28372

$`40586`
 [1]     1     2     7    13    14    26    91   182   223   446
[11]  1561  2899  3122  5798 20293 40586

$`51390`
 [1]     1     2     3     5     6     9    10    15    18    30
[11]    45    90   571  1142  1713  2855  3426  5139  5710  8565
[21] 10278 17130 25695 51390

$`98248`
[1]     1     2     4     8 12281 24562 49124 98248


## Quickly generate large primes over small interval
options(scipen = 50)
system.time(myPs <- primeSieve(10^13+10^3, 10^13))
   user  system elapsed 
  0.016   0.000   0.016
  
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

Contact
----
I welcome any and all feedback. If you would like to report a bug, have a question, or have suggestions for possible improvements, please contact me here: jwood000@gmail.com
