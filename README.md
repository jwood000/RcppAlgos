# RcppAlgos

Overview
---------
A collection of high performance functions implemented in C++ with Rcpp for solving problems in combinatorics and computational mathematics. Featured functions:

* `primeSieve` - Generates all primes less than a billion in just over 1 second
* `primeCount` -  Counts the number of primes below a trillion in under 0.5 seconds.
* `comboGeneral`/`permuteGeneral` - Generate all combinations/permutations of a vector (including [multisets](https://en.wikipedia.org/wiki/Multiset)) meeting specific criteria.
    - Produce results in parallel using the `Parallel` argument (**development version only**). You can also apply each of the five compiled functions given by the argument `constraintFun` in parallel as well. E.g. Obtaining the row sums of all combinations:
        - `comboGeneral(20, 10, constraintFun = "sum", Parallel = TRUE)`
    - Alternatively, the arguments `lower` and `upper` make it possible to generate combinations/permutations in chunks allowing for parallelization via the package `parallel`. This is convenient when you want to apply a custom function to the output in parallel as well (see this [stackoverflow post](https://stackoverflow.com/a/51595866/4408538) for a use case).
    - GMP support allows for exploration of combinations/permutations of vectors with many elements.
* `comboSample`/`permuteSample` - Easily generate random samples of combinations/permutations in parallel.
    - You can pass a vector of specific indices or rely on the internal sampling functions. We call `sample` when the total number of results is small and for larger cases, the sampling is done in a very similar fashion to `urand.bigz` from the `gmp` package.

The `primeSieve` function and the `primeCount` function are both based off of the excellent work by [Kim Walisch](https://github.com/kimwalisch). The respective repos can be found here: [kimwalisch/primesieve](https://github.com/kimwalisch/primesieve); [kimwalisch/primecount](https://github.com/kimwalisch/primecount)

Additionally, many of the sieving functions make use of the fast integer division library [libdivide](https://github.com/ridiculousfish/libdivide) by [ridiculousfish](https://github.com/ridiculousfish).

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
Easily executed with a very simple interface. Output is in [lexicographical order](https://en.wikipedia.org/wiki/Lexicographical_order).
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
 [2,]    1    2    4
 [3,]    1    3    2
 [4,]    1    3    4
  .      .    .    .
  .      .    .    .
[21,]    4    2    1
[22,]    4    2    3
[23,]    4    3    1
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
  0.124   0.058   0.182 

nrow(comboGeneral(25,13))
[1] 5200300
```

### Using Constraints
Oftentimes, one needs to find combinations/permutations that meet certain requirements.

``` r
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
 
tail(p)
        [,1] [,2] [,3] [,4] [,5] [,6]
[3955,]  354  287  149  187  199 1176
[3956,]  354  287  149  199  187 1176
[3957,]  354  287  187  149  199 1176
[3958,]  354  287  187  199  149 1176
[3959,]  354  287  199  149  187 1176
[3960,]  354  287  199  187  149 1176


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
```

### Working with Multisets
Sometimes, the standard combination/permutation functions don't quite get us to our desired goals. For
example, one may need all permutations of a vector with some of the elements repeated a specific
amount of times (i.e. a multiset). Consider the following vector `a <- c(1,1,1,1,2,2,2,7,7,7,7,7)` and one
would like to find permutations of `a` of length 6. Using traditional methods, we would need to generate all
permutations, then eliminate duplicate values. Even considering that `permuteGeneral` is very efficient,
this approach is clunky and not as fast as it could be. Observe:

``` r
getPermsWithSpecificRepetition <- function(z, n) {
    b <- permuteGeneral(z, n)
    myDupes <- duplicated(apply(b, 1, paste, collapse=""))
    b[!myDupes, ]
}

system.time(test <- getPermsWithSpecificRepetition(a, 6))
   user  system elapsed 
  4.300   0.028   4.331
```


#### Enter _freqs_
Situations like this call for the use of the `freqs` argument. Simply, enter the number
of times each unique element is repeated and Voila!

``` r
system.time(test2 <- permuteGeneral(unique(a), 6, freqs = rle(a)$lengths))
   user  system elapsed 
      0       0       0
      
identical(test, test2)
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

## or combinations of a certain length
comboGeneral(3, 2, freqs = c(1,2,2), constraintFun = "prod")
     [,1] [,2] [,3]
[1,]    1    2    2
[2,]    1    3    3
[3,]    2    2    4
[4,]    2    3    6
[5,]    3    3    9
```

##### All Combinatorial Functions Work with Factors

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

### Parallel Computing
Using the parameter `Parallel` (**development version only**), we can easily generate combinations/permutations with great efficiency.

```r
## RcppAlgos uses the number of cores available minus one
parallel::detectCores()
[1] 8

identical(comboGeneral(20, 10, freqs = rep(1:4, 5)),
          comboGeneral(20, 10, freqs = rep(1:4, 5), Parallel = TRUE))
[1] TRUE

## Using 7 cores
library(microbenchmark)
microbenchmark(serial = comboGeneral(20, 10, freqs = rep(1:4, 5)),
             parallel = comboGeneral(20, 10, freqs = rep(1:4, 5), Parallel = TRUE))
Unit: milliseconds
     expr       min        lq      mean    median        uq      max neval
   serial 216.08683 235.09945 241.97200 240.59142 247.37669 299.7269   100
 parallel  60.10095  64.24945  76.92387  77.25798  85.44578 118.6835   100
```

And applying any of the constraint functions in parallel is highly efficient as well. Consider obtaining the row sums of all combinations:

```r
## base R using combn and FUN
combnSum <- combn(20, 10, sum)
algosSum <- comboGeneral(20, 10, constraintFun = "sum")

identical(as.integer(combnSum), algosSum[,11])
[1] TRUE

## Using parallel
paralSum <- comboGeneral(20, 10, constraintFun = "sum", Parallel = TRUE)
identical(paralSum, algosSum)
[1] TRUE

microbenchmark(serial = comboGeneral(20, 10, constraintFun = "sum"),
             parallel = comboGeneral(20, 10, constraintFun = "sum", Parallel = TRUE),
             combnSum = combn(20, 10, sum))
Unit: milliseconds
     expr        min         lq       mean     median         uq        max neval
   serial   3.305383   3.943742   4.340827   3.998379   4.175258   7.638484   100
 parallel   1.042070   1.391787   1.484434   1.433409   1.479612   3.277695   100
 combnSum 203.691824 215.036510 220.745091 217.890390 222.273996 304.720899   100
```
### Faster than `rowSums` and `rowMeans`
In fact, finding row sums or row means is even faster than simply applying the highly efficient `rowSums`/`rowMeans` after the combinations have already been generated:
```r
## Pre-generate combinations
combs <- comboGeneral(25, 10)

## Testing rowSums alone against generating combinations as well as summing
microbenchmark(serial = comboGeneral(25, 10, constraintFun = "sum"),
             parallel = comboGeneral(25, 10, constraintFun = "sum", Parallel = TRUE),
              rowsums = rowSums(combs))
Unit: milliseconds
     expr       min        lq      mean    median        uq      max neval
   serial 112.79266 117.04039 126.27396 119.29226 122.72389 200.1609   100
 parallel  39.10258  41.95387  51.79792  46.31817  49.83576 115.4092   100
  rowsums 103.22926 104.30309 109.73375 105.28372 111.05051 183.4639   100

all.equal(rowSums(combs), 
          comboGeneral(25, 10, 
                       constraintFun = "sum",
                       Parallel = TRUE)[,11])
[1] TRUE

## Testing rowMeans alone against generating combinations as well as obtain row means
microbenchmark(serial = comboGeneral(25, 10, constraintFun = "mean"),
             parallel = comboGeneral(25, 10, constraintFun = "mean", Parallel = TRUE),
             rowmeans = rowMeans(combs))
Unit: milliseconds
     expr       min        lq      mean    median       uq      max neval
   serial 173.58781 183.90069 211.03830 196.17403 241.7987 307.3608   100
 parallel  53.10161  57.80737  77.03166  62.58454 103.2152 228.0634   100
 rowmeans 108.65154 111.16465 117.70982 114.76217 123.3963 136.3974   100
 
all.equal(rowMeans(combs), 
          comboGeneral(25, 10, 
                       constraintFun = "mean",
                       Parallel = TRUE)[,11])
[1] TRUE
```
We are doing double the work nearly twice as fast in both cases!!

### Using arguments `lower` and `upper`
There are arguments `lower` and `upper` that can be utilized to generate chunks of combinations/permutations without having to generate all of them followed by subsetting.  As the output is in lexicographical order, these arguments specify where to start and stop generating. For example, `comboGeneral(5, 3)` outputs 10 combinations of the vector `1:5` chosen 3 at a time. We can set `lower` to 5 in order to start generation from the 5<sup>th</sup> lexicographical combination. Similarly, we can set `upper` to 4 in order only generate the first 4 combinations. We can also use them together to produce only a certain chunk of combinations. For example, setting `lower` to 4 and `upper` to 6 only produces the 4<sup>th</sup>, 5<sup>th</sup>, and 6<sup>th</sup> lexicographical combinations. Observe:

``` r
comboGeneral(5, 3, lower = 4, upper = 6)
     [,1] [,2] [,3]
[1,]    1    3    4
[2,]    1    3    5
[3,]    1    4    5

## is equivalent to the following:
comboGeneral(5, 3)[4:6, ]
     [,1] [,2] [,3]
[1,]    1    3    4
[2,]    1    3    5
[3,]    1    4    5
```
In addition to being useful by avoiding the unnecessary overhead of generating all combination/permutations followed by subsetting just to see a few specific results, `lower` and `upper` can be utilized to generate large number of combinations/permutations in parallel. Observe:

``` r
## Over 3 billion results
comboCount(35, 15)
[1] 3247943160

## 10086780 evenly divides 3247943160

system.time(lapply(seq(1, 3247943160, 10086780), function(x) {
        temp <- comboGeneral(35, 15, lower = x, upper = x + 10086779)
        ## do something
        x
}))
   user  system elapsed 
 99.495  49.993 149.467

## Enter parallel

library(parallel)
system.time(mclapply(seq(1, 3247943160, 10086780), function(x) {
     temp <- comboGeneral(35, 15, lower = x, upper = x + 10086779)
     ## do something
     x
}, mc.cores = 6))
   user  system elapsed 
125.361  85.916  35.641
```
These arguments are also useful when one needs to explore combinations/permutations of really large vectors:
```r
set.seed(222)
myVec <- rnorm(1000)

## HUGE number of combinations
comboCount(myVec, 50, repetition = TRUE)
Big Integer ('bigz') :
[1] 109740941767310814894854141592555528130828577427079559745647393417766593803205094888320

## Let's look at one hundred thousand combinations in the range (1e15 + 1, 1e15 + 1e5)
system.time(b <- comboGeneral(myVec, 50, TRUE, 
                              lower = 1e15 + 1,
                              upper = 1e15 + 1e5))
   user  system elapsed 
  0.008   0.001   0.010
  
b[1:5, 45:50]
          [,1]      [,2]      [,3]     [,4]      [,5]       [,6]
[1,] 0.5454861 0.4787456 0.7797122 2.004614 -1.257629 -0.7740501
[2,] 0.5454861 0.4787456 0.7797122 2.004614 -1.257629  0.1224679
[3,] 0.5454861 0.4787456 0.7797122 2.004614 -1.257629 -0.2033493
[4,] 0.5454861 0.4787456 0.7797122 2.004614 -1.257629  1.5511027
[5,] 0.5454861 0.4787456 0.7797122 2.004614 -1.257629  1.0792094
```
### Sampling
We can also produce random samples of combinations/permutations with `comboSample` and `permuteSample`. This is really useful when we need a reproducible set of random combinations/permutations. Many of the traditional ways of doing this involved relying on heavy use of `sample` and hoping that we don't generate duplicate results. Both functions have a similar interface to their respective `General` functions. Observe:

``` r
comboSample(10, 8, TRUE, n = 5, seed = 84)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
[1,]    3    3    3    6    6   10   10   10
[2,]    1    3    3    4    4    7    9   10
[3,]    3    7    7    7    9   10   10   10
[4,]    3    3    3    9   10   10   10   10
[5,]    1    2    2    3    3    4    4    7

## We can also use sampleVec to generate specific results
## E.g. the below generates the 1st, 5th, 25th, 125th, and
## 625th lexicographical combinations
comboSample(10, 8, TRUE, sampleVec = c(1, 5, 25, 125, 625))
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
[1,]    1    1    1    1    1    1    1    1
[2,]    1    1    1    1    1    1    1    5
[3,]    1    1    1    1    1    1    3    8
[4,]    1    1    1    1    1    3    6    9
[5,]    1    1    1    1    5    6   10   10

## Is the same as:
comboGeneral(10, 8, TRUE)[5^(0:4), ]
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
[1,]    1    1    1    1    1    1    1    1
[2,]    1    1    1    1    1    1    1    5
[3,]    1    1    1    1    1    1    3    8
[4,]    1    1    1    1    1    3    6    9
[5,]    1    1    1    1    5    6   10   10
```
Just like the `General` counterparts (i.e. `combo/permuteGeneral`), we can easily explore combinations/permutations of large vectors where the total number of results is enormous in parallel (**development version only**).
```r
permuteSample(500, 10, TRUE, n = 5, seed = 123, Parallel = TRUE)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
[1,]   55  435  274  324  200  152    6  313  121   377
[2,]  196  166  331  154  443  329  155  233  354   442
[3,]  235  325   94   27  370  117  302   86  229   126
[4,]  284  104  464  104  207  127  117    9  390   414
[5,]  456   76  381  456  219   23  376  187   11   123

permuteSample(factor(state.abb), 15, n = 3, seed = 50, Parallel = TRUE)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15]
[1,] ME   FL   DE   OK   ND   CA   PA   AL   ID   MO    NM    HI    KY    MT    NJ   
[2,] AZ   CA   AL   CT   ME   SD   ID   SC   OK   NH    HI    TN    ND    IA    MT   
[3,] MD   MO   NC   MT   NH   AL   VA   MA   VT   WV    NJ    NE    MN    MS    MI   
50 Levels: AK AL AR AZ CA CO CT DE FL GA HI IA ID IL IN KS KY LA MA MD ME MI MN ... WY

permuteCount(factor(state.abb), 15)
Big Integer ('bigz') :
[1] 2943352142120754524160000
```
### User Defined Functions
You can also pass user defined functions by utilizing the argument `FUN`. This feature's main purpose is for convenience, however it is somewhat more efficient than generating all combinations/permutations and then using a function from the `apply` family (N.B. the argument `Parallel` has no effect when `FUN` is employed).

```r
funCustomComb = function(n, r) {
    combs = comboGeneral(n, r)
    lapply(1:nrow(combs), function(x) cumprod(combs[x,]))
}

identical(funCustomComb(15, 8), comboGeneral(15, 8, FUN = cumprod))
[1] TRUE

library(microbenchmark)
microbenchmark(f1 = funCustomComb(15, 8),
                f2 = comboGeneral(15, 8, FUN = cumprod), unit = "relative")
unit: relative
 expr      min       lq     mean   median       uq      max neval
   f1 6.946481 6.891553 6.334866 6.821221 6.934111 2.686777   100
   f2 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000   100
   
comboGeneral(15, 8, FUN = cumprod, upper = 3)
[[1]]
[1]     1     2     6    24   120   720  5040 40320

[[2]]
[1]     1     2     6    24   120   720  5040 45360

[[3]]
[1]     1     2     6    24   120   720  5040 50400

## Works with the sampling functions as well
permuteSample(5000, 1000, n = 3, seed = 101, FUN = sd)
[[1]]
[1] 1431.949

[[2]]
[1] 1446.859

[[3]]
[1] 1449.272
```

Mathematical Computation
-----
`RcppAlgos` comes equipped with several functions for quickly generating essential components
for problems common in computational mathematics.

The following sieving functions (`primeFactorizeSieve`, `divisorsSieve`, `numDivisorSieve`, & `eulerPhiSieve`) are very useful and flexible. Generate components up to a number or between two bounds.

``` r
## get the number of divisors for every number from 1 to n
numDivisorSieve(20)
 [1] 1 2 2 3 2 4 2 4 3 4 2 6 2 4 4 5 2 6 2 6

## If you want the complete factorization from 1 to n, use divisorsList
system.time(allFacs <- divisorsSieve(10^5, namedList = TRUE))
   user  system elapsed 
  0.073   0.004   0.077

allFacs[c(4339, 15613, 22080)]
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


## Between two bounds
primeFactorizeSieve(10^12, 10^12 + 5)
[[1]]
 [1] 2 2 2 2 2 2 2 2 2 2 2 2 5 5 5 5 5 5 5 5 5 5 5 5

[[2]]
[1]       73      137 99990001

[[3]]
[1]            2            3 166666666667

[[4]]
[1]      61   14221 1152763

[[5]]
[1]      2      2     17    149    197 501001

[[6]]
[1]           3           5 66666666667


## Creating a named object
eulerPhiSieve(20, namedVector = TRUE)
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
 1  1  2  2  4  2  6  4  6  4 10  4 12  6  8  8 16  6 18  8
```

### Vectorized Functions
There are three very fast vectorized functions for general factoring (e.g. all divisors of number), primality testing, as well as prime factoring (`divisorsRcpp`, `isPrimeRcpp`, `primeFactorize`)

``` r
## get result for individual numbers
primeFactorize(123456789)
[1]    3    3 3607 3803


## or for an entire vector
set.seed(100)
myVec <- sample(-100000000:100000000, 5)
divisorsRcpp(myVec, namedList = TRUE)
$`-38446778`
[1] -38446778 -19223389        -2        -1         1
[6]         2  19223389  38446778

$`-48465500`
 [1] -48465500 -24232750 -12116375  -9693100  -4846550
 [6]  -2423275  -1938620   -969310   -484655   -387724
[11]   -193862    -96931      -500      -250      -125
[16]      -100       -50       -25       -20       -10
[21]        -5        -4        -2        -1         1
[26]         2         4         5        10        20
[31]        25        50       100       125       250
[36]       500     96931    193862    387724    484655
[41]    969310   1938620   2423275   4846550   9693100
[46]  12116375  24232750  48465500

$`10464487`
[1]        1       11      317     3001     3487    33011
[7]   951317 10464487

$`-88723370`
 [1] -88723370 -44361685 -17744674  -8872337       -10
 [6]        -5        -2        -1         1         2
[11]         5        10   8872337  17744674  44361685
[16]  88723370

$`-6290143`
[1] -6290143       -1        1  6290143


## Creating a named object
isPrimeRcpp(995:1000, namedVector = TRUE)
  995   996   997   998   999  1000 
FALSE FALSE  TRUE FALSE FALSE FALSE

```


### primeSieve & primeCount
Both of these functions are based on the excellent algorithms developed by [Kim Walisch](https://github.com/kimwalisch).

``` r
## Quickly generate large primes over small interval
options(scipen = 50)
system.time(myPs <- primeSieve(10^13+10^3, 10^13))
   user  system elapsed 
  0.035   0.002   0.036
  
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

## primes under a billion!!!
system.time(primeSieve(10^9))
   user  system elapsed 
  1.289   0.091   1.382
  
  
## Enumerate the number of primes below trillion
system.time(underOneTrillion <- primeCount(10^12))
   user  system elapsed 
  0.484   0.000   0.485
  
underOneTrillion
[1] 37607912018


## Enumerate the number of primes below a billion in 2 milliseconds
library(microbenchmark)
microbenchmark(primeCount(10^9))
Unit: milliseconds
             expr      min      lq    mean  median       uq      max neval
 primeCount(10^9) 1.998919 2.00191 2.08563 2.00397 2.082288 3.043241   100
 
 
primeCount(10^9)
[1] 50847534
```

Contact
----
I welcome any and all feedback. If you would like to report a bug, have a question, or have suggestions for possible improvements, please contact me here: jwood000@gmail.com
