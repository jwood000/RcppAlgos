# High Performance Benchmarks

This document serves as an overview for measuring the performance of
`RcppAlgos` against other tools for generating combinations,
permutations, and partitions. This stackoverflow post: [How to generate
permutations or combinations of object in
R?](https://stackoverflow.com/a/47983855/4408538) has some benchmarks.
You will note that the examples in that post are relatively small. The
benchmarks below will focus on larger examples where performance really
matters and for this reason we only consider the packages
[arrangements](https://cran.r-project.org/package=arrangements),
[partitions](https://cran.r-project.org/package=partitions), and
[RcppAlgos](https://cran.r-project.org/package=RcppAlgos).

## Setup Information

For the benchmarks below, we used a `2025 Macbook Air Apple M4 24 GB`
machine.

``` r

library(RcppAlgos)
library(partitions)
library(arrangements)
#> 
#> Attaching package: 'arrangements'
#> The following object is masked from 'package:partitions':
#> 
#>     compositions
library(microbenchmark)

options(digits = 4)
options(width = 90)

pertinent_output <- capture.output(sessionInfo())
cat(paste(pertinent_output[1:3], collapse = "\n"))
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Sequoia 15.7.4

pkgs <- c("RcppAlgos", "arrangements", "partitions", "microbenchmark")
sapply(pkgs, packageVersion, simplify = FALSE)
#> $RcppAlgos
#> [1] '2.10.3'
#> 
#> $arrangements
#> [1] '1.1.10'
#> 
#> $partitions
#> [1] '1.10.9'
#> 
#> $microbenchmark
#> [1] '1.5.0'

numThreads <- min(as.integer(RcppAlgos::stdThreadMax() / 2), 6)
numThreads
#> [1] 5
```

## Combinations

### Combinations - Distinct

``` r

set.seed(13)
v1 <- sort(sample(100, 30))
m <- 21
t1 <- comboGeneral(v1, m, Parallel = T)
t2 <- combinations(v1, m)
stopifnot(identical(t1, t2))
dim(t1)
#> [1] 14307150       21
rm(t1, t2)
invisible(gc())
microbenchmark(cbRcppAlgosPar = comboGeneral(v1, m, nThreads = numThreads),
               cbRcppAlgosSer = comboGeneral(v1, m),
               cbArrangements = combinations(v1, m),
               times = 15, unit = "relative")
#> Warning in microbenchmark(cbRcppAlgosPar = comboGeneral(v1, m, nThreads = numThreads), :
#> less accurate nanosecond times to avoid potential integer overflows
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.000    15
#>  cbRcppAlgosSer 3.546 3.018 2.864  2.786 2.701 2.635    15
#>  cbArrangements 2.720 2.285 2.149  2.099 1.995 1.954    15
```

### Combinations - Repetition

``` r

v2 <- v1[1:10]
m <- 20
t1 <- comboGeneral(v2, m, repetition = TRUE, nThreads = numThreads)
t2 <- combinations(v2, m, replace = TRUE)
stopifnot(identical(t1, t2))
dim(t1)
#> [1] 10015005       20
rm(t1, t2)
invisible(gc())
microbenchmark(cbRcppAlgosPar = comboGeneral(v2, m, TRUE, nThreads = numThreads),
               cbRcppAlgosSer = comboGeneral(v2, m, TRUE),
               cbArrangements = combinations(v2, m, replace = TRUE),
               times = 15, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.000    15
#>  cbRcppAlgosSer 3.053 2.984 2.799  2.869 2.641 2.573    15
#>  cbArrangements 2.096 2.112 1.972  2.021 1.868 1.811    15
```

### Combinations - Multisets

``` r

myFreqs <- c(2, 4, 4, 5, 3, 2, 2, 2, 3, 4, 1, 4, 2, 5)
v3 <- as.integer(c(1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610))
t1 <- comboGeneral(v3, 20, freqs = myFreqs, nThreads = numThreads)
t2 <- combinations(freq = myFreqs, k = 20, x = v3)
stopifnot(identical(t1, t2))
dim(t1)
#> [1] 14594082       20
rm(t1, t2)
invisible(gc())
microbenchmark(cbRcppAlgosPar = comboGeneral(v3, 20, freqs = myFreqs, nThreads = numThreads),
               cbRcppAlgosSer = comboGeneral(v3, 20, freqs = myFreqs),
               cbArrangements = combinations(freq = myFreqs, k = 20, x = v3),
               times = 10, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.000    10
#>  cbRcppAlgosSer 2.978 3.160 3.178  3.045 2.839 4.554    10
#>  cbArrangements 5.130 5.107 4.910  5.000 4.723 4.638    10
```

## Permutations

### Permutations - Distinct

``` r

v4 <- as.integer(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59))
t1 <- permuteGeneral(v4, 6, nThreads = numThreads)
t2 <- permutations(v4, 6)
stopifnot(identical(t1, t2))
dim(t1)
#> [1] 8910720       6
rm(t1, t2)
invisible(gc())
microbenchmark(cbRcppAlgosPar = permuteGeneral(v4, 6, nThreads = numThreads),
               cbRcppAlgosSer = permuteGeneral(v4, 6),
               cbArrangements = permutations(v4, 6),
               times = 15, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.000    15
#>  cbRcppAlgosSer 1.732 1.701 1.535  1.722 1.570 1.143    15
#>  cbArrangements 2.335 2.517 2.371  2.485 2.605 2.118    15


## Indexing permutation example with the partitions package
t1 <- permuteGeneral(11, nThreads = 4)
t2 <- permutations(11)
t3 <- perms(11)

dim(t1)
#> [1] 39916800       11

stopifnot(identical(t1, t2), identical(t1, t(as.matrix(t3))))
rm(t1, t2, t3)
invisible(gc())

microbenchmark(cbRcppAlgosPar = permuteGeneral(11, nThreads = 4),
               cbRcppAlgosSer = permuteGeneral(11),
               cbArrangements = permutations(11),
               cbPartitions   = perms(11),
               times = 5, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.000     5
#>  cbRcppAlgosSer 2.714 2.529 2.542  2.430 2.248 2.795     5
#>  cbArrangements 3.970 3.900 3.758  3.776 3.412 3.801     5
#>    cbPartitions 8.680 8.363 8.927  9.611 8.607 9.281     5
```

### Permutations - Repetition

``` r

v5 <- v3[1:5]
t1 <- permuteGeneral(v5, 10, repetition = TRUE, nThreads = numThreads)
t2 <- permutations(v5, 10, replace = TRUE)
stopifnot(identical(t1, t2))
dim(t1)
#> [1] 9765625      10
rm(t1, t2)
invisible(gc())
microbenchmark(cbRcppAlgosPar = permuteGeneral(v5, 10, TRUE, nThreads = numThreads),
               cbRcppAlgosSer = permuteGeneral(v5, 10, TRUE),
               cbArrangements = permutations(x = v5, k = 10, replace = TRUE),
               times = 10, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq    max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.0000    10
#>  cbRcppAlgosSer 2.471 2.475 2.098  2.646 2.506 0.9127    10
#>  cbArrangements 2.508 2.708 2.331  2.711 2.711 1.4020    10
```

### Permutations - Multisets

``` r

v6 <- sort(runif(12))
t1 <- permuteGeneral(v6, 7, freqs = rep(1:3, 4), nThreads = numThreads)
t2 <- permutations(freq = rep(1:3, 4), k = 7, x = v6)
stopifnot(identical(t1, t2))
dim(t1)
#> [1] 19520760        7
rm(t1, t2)
invisible(gc())
microbenchmark(cbRcppAlgosPar = permuteGeneral(v6, 7, freqs = rep(1:3, 4), nThreads = numThreads),
               cbRcppAlgosSer = permuteGeneral(v6, 7, freqs = rep(1:3, 4)),
               cbArrangements = permutations(freq = rep(1:3, 4), k = 7, x = v6),
               times = 10, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.000    10
#>  cbRcppAlgosSer 3.430 3.334 4.241  4.078 4.901 4.911    10
#>  cbArrangements 3.666 3.549 3.980  3.791 4.376 4.471    10
```

## Partitions

### Partitions - Distinct

#### All Distinct Partitions

``` r

t1 <- comboGeneral(0:140, freqs=c(140, rep(1, 140)),
                   constraintFun = "sum", comparisonFun = "==",
                   limitConstraints = 140)
t2 <- partitions(140, distinct = TRUE)
t3 <- diffparts(140)

# Each package has different output formats... we only examine dimensions
# and that each result is a partition of 140
stopifnot(identical(dim(t1), dim(t2)), identical(dim(t1), dim(t(t3))),
                    all(rowSums(t1) == 140), all(rowSums(t2) == 140),
                    all(colSums(t3) == 140))
dim(t1)
#> [1] 9617150      16
rm(t1, t2, t3)
invisible(gc())
microbenchmark(cbRcppAlgosPar = partitionsGeneral(0:140, freqs=c(140, rep(1, 140)), nThreads = numThreads),
               cbRcppAlgosSer = partitionsGeneral(0:140, freqs=c(140, rep(1, 140))),
               cbArrangements = partitions(140, distinct = TRUE),
               cbPartitions   = diffparts(140),
               times = 10, unit = "relative")
#> Unit: relative
#>            expr    min     lq   mean median     uq    max neval
#>  cbRcppAlgosPar  1.000  1.000  1.000  1.000  1.000  1.000    10
#>  cbRcppAlgosSer  3.264  2.742  2.638  2.471  2.491  2.442    10
#>  cbArrangements  2.934  2.519  2.382  2.243  2.197  2.186    10
#>    cbPartitions 17.438 14.845 13.948 13.329 12.822 12.759    10
```

#### Restricted Distinct Partitions

``` r

t1 <- comboGeneral(160, 10,
                   constraintFun = "sum", comparisonFun = "==",
                   limitConstraints = 160)
t2 <- partitions(160, 10, distinct = TRUE)
stopifnot(identical(t1, t2))
dim(t1)
#> [1] 8942920      10
rm(t1, t2)
invisible(gc())
microbenchmark(cbRcppAlgosPar = partitionsGeneral(160, 10, nThreads = numThreads),
               cbRcppAlgosSer = partitionsGeneral(160, 10),
               cbArrangements = partitions(160, 10, distinct = TRUE),
               times = 10, unit = "relative")
#> Unit: relative
#>            expr   min    lq mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.00  1.000 1.000 1.000    10
#>  cbRcppAlgosSer 3.065 3.082 3.16  3.190 3.349 2.996    10
#>  cbArrangements 3.776 3.750 3.74  3.692 3.768 3.499    10
```

### Partitions - Repetition

#### All Partitions

``` r

t1 <- comboGeneral(0:65, repetition = TRUE, constraintFun = "sum",
                   comparisonFun = "==", limitConstraints = 65)
t2 <- partitions(65)
t3 <- parts(65)

# Each package has different output formats... we only examine dimensions
# and that each result is a partition of 65
stopifnot(identical(dim(t1), dim(t2)), identical(dim(t1), dim(t(t3))),
          all(rowSums(t1) == 65), all(rowSums(t2) == 65),
          all(colSums(t3) == 65))
dim(t1)
#> [1] 2012558      65
rm(t1, t2, t3)
invisible(gc())
microbenchmark(cbRcppAlgosPar = partitionsGeneral(0:65, repetition = TRUE,
                                                  nThreads = numThreads),
               cbRcppAlgosSer = partitionsGeneral(0:65, repetition = TRUE),
               cbArrangements = partitions(65),
               cbPartitions   = parts(65),
               times = 20, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.000    20
#>  cbRcppAlgosSer 2.768 2.765 2.185  2.623 2.114 1.664    20
#>  cbArrangements 2.788 2.779 2.143  2.649 2.149 1.312    20
#>    cbPartitions 8.427 8.466 6.428  8.218 6.218 3.945    20
```

#### Restricted Partitions

``` r

t1 <- comboGeneral(100, 15, TRUE, constraintFun = "sum",
                   comparisonFun = "==", limitConstraints = 100)
t2 <- partitions(100, 15)
stopifnot(identical(t1, t2))
dim(t1)
#> [1] 9921212      15
rm(t1, t2)

# This takes a really long time... not because of restrictedparts,
# but because apply is not that fast. This transformation is
# needed for proper comparisons. As a result, we will compare
# a smaller example
# t3 <- t(apply(as.matrix(restrictedparts(100, 15, include.zero = F)), 2, sort))
t3 <- t(apply(as.matrix(restrictedparts(50, 15, include.zero = F)), 2, sort))
stopifnot(identical(partitions(50, 15), t3))
rm(t3)
invisible(gc())
microbenchmark(cbRcppAlgosPar = partitionsGeneral(100, 15, TRUE,
                                                  nThreads = numThreads),
               cbRcppAlgosSer = partitionsGeneral(100, 15, TRUE),
               cbArrangements = partitions(100, 15),
               cbPartitions   = restrictedparts(100, 15,
                                                include.zero = FALSE),
               times = 10, unit = "relative")
#> Unit: relative
#>            expr    min     lq   mean median     uq   max neval
#>  cbRcppAlgosPar  1.000  1.000  1.000  1.000  1.000 1.000    10
#>  cbRcppAlgosSer  3.326  3.107  2.732  2.812  2.435 2.307    10
#>  cbArrangements  3.473  3.276  2.890  3.111  2.528 2.384    10
#>    cbPartitions 14.466 13.713 11.828 12.445 10.159 9.943    10
```

### Partitions - Multisets

Currently, `RcppAlgos` is the only package capable of efficiently
generating partitions of multisets. Therefore, we will only time
`RcppAlgos` and use this as a reference for future improvements.

``` r

t1 <- comboGeneral(120, 10, freqs = rep(1:8, 15),
                   constraintFun = "sum", comparisonFun = "==",
                   limitConstraints = 120)
dim(t1)
#> [1] 7340225      10
stopifnot(all(rowSums(t1) == 120))
microbenchmark(cbRcppAlgos = partitionsGeneral(120, 10, freqs = rep(1:8, 15)),
               times = 10)
#> Unit: milliseconds
#>         expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgos 88.17 88.87 92.39  89.26 89.96 115.6    10
```

## Compositions

### Compositions - Repetition

#### All Compositions (Small case)

``` r

t1 <- compositionsGeneral(0:15, repetition = TRUE)
t2 <- arrangements::compositions(15)
t3 <- partitions::compositions(15)

# Each package has different output formats... we only examine dimensions
# and that each result is a partition of 15
stopifnot(identical(dim(t1), dim(t2)), identical(dim(t1), dim(t(t3))),
          all(rowSums(t1) == 15), all(rowSums(t2) == 15),
          all(colSums(t3) == 15))
dim(t1)
#> [1] 16384    15
rm(t1, t2, t3)
invisible(gc())
microbenchmark(cbRcppAlgosSer = compositionsGeneral(0:15, repetition = TRUE),
               cbArrangements = arrangements::compositions(15),
               cbPartitions   = partitions::compositions(15),
               times = 20, unit = "relative")
#> Unit: relative
#>            expr     min      lq   mean  median      uq     max neval
#>  cbRcppAlgosSer   1.008   1.002   1.01   1.017   1.014   1.035    20
#>  cbArrangements   1.000   1.000   1.00   1.000   1.000   1.000    20
#>    cbPartitions 112.405 122.651 151.68 148.691 182.256 190.316    20
```

For the next two examples, we will exclude the `partitions` package for
efficiency reasons.

#### All Compositions (Larger case)

``` r

t1 <- compositionsGeneral(0:23, repetition = TRUE)
t2 <- arrangements::compositions(23)

# Each package has different output formats... we only examine dimensions
# and that each result is a partition of 23
stopifnot(identical(dim(t1), dim(t2)), all(rowSums(t1) == 23),
          all(rowSums(t2) == 23))
dim(t1)
#> [1] 4194304      23
rm(t1, t2)
invisible(gc())
microbenchmark(cbRcppAlgosPar = compositionsGeneral(0:23, repetition = TRUE,
                                                    nThreads = numThreads),
               cbRcppAlgosSer = compositionsGeneral(0:23, repetition = TRUE),
               cbArrangements = arrangements::compositions(23),
               times = 20, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.000    20
#>  cbRcppAlgosSer 3.684 3.658 3.713  3.734 3.776 3.730    20
#>  cbArrangements 3.618 3.573 3.575  3.586 3.575 3.549    20
```

#### Compositions of Specific Length

``` r

t1 <- compositionsGeneral(30, 10, repetition = TRUE)
t2 <- arrangements::compositions(30, 10)

stopifnot(identical(t1, t2), all(rowSums(t1) == 30))
dim(t1)
#> [1] 10015005       10
rm(t1, t2)
invisible(gc())
microbenchmark(cbRcppAlgosPar = compositionsGeneral(30, 10, repetition = TRUE,
                                                    nThreads = numThreads),
               cbRcppAlgosSer = compositionsGeneral(30, 10, repetition = TRUE),
               cbArrangements = arrangements::compositions(30, 10),
               times = 20, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosPar 1.000 1.000 1.000  1.000 1.000 1.000    20
#>  cbRcppAlgosSer 2.892 2.836 2.642  2.818 2.602 1.486    20
#>  cbArrangements 3.264 3.101 2.886  3.080 2.833 1.685    20
```

### Specialized Composition Benchmarks

Similar to partitions of multisets, several composition variants
supported in `RcppAlgos` (e.g., distinct parts, capped cases, etc.) are
not currently available in other R packages or combinatorics libraries.
In such cases, benchmarks will only report timings for `RcppAlgos`,
which can serve as a reference point for future implementations and
improvements.

#### Compositions with Specific `target`

``` r

t1 <- compositionsGeneral(10, 10, repetition = TRUE, target = 30)
dim(t1)
#> [1] 9091270      10

stopifnot(all(rowSums(t1) == 30))
microbenchmark(
    cbRcppAlgosSer = compositionsGeneral(10, 10, repetition = TRUE, target = 30),
    cbRcppAlgosPar = compositionsGeneral(10, 10, repetition = TRUE,
                                         target = 30, nThreads = numThreads),
    times = 10,
    check = "identical"
)
#> Unit: milliseconds
#>            expr   min    lq  mean median    uq   max neval
#>  cbRcppAlgosSer 57.01 58.14 59.83  59.34 62.12 62.44    10
#>  cbRcppAlgosPar 17.52 19.75 23.58  22.45 23.89 38.51    10
```

#### Compositions with Distinct Parts

``` r

t1 <- compositionsGeneral(50, 8)
dim(t1)
#> [1] 4677120       8

stopifnot(all(rowSums(t1) == 50))
microbenchmark(
    cbRcppAlgosSer = compositionsGeneral(50, 8),
    cbRcppAlgosPar = compositionsGeneral(50, 8, nThreads = numThreads),
    times = 10,
    check = "identical"
)
#> Unit: milliseconds
#>            expr    min     lq   mean median     uq    max neval
#>  cbRcppAlgosSer 214.55 215.45 216.71 216.19 217.17 223.18    10
#>  cbRcppAlgosPar  54.08  55.45  57.64  56.84  59.71  62.19    10
```

#### Compositions with Distinct Parts & Specific `target`

``` r

t1 <- compositionsGeneral(30, 7, target = 60)
dim(t1)
#> [1] 11773440        7

stopifnot(all(rowSums(t1) == 60))
microbenchmark(
    cbRcppAlgosSer = compositionsGeneral(30, 7, target = 60),
    cbRcppAlgosPar = compositionsGeneral(30, 7, target = 60, nThreads = numThreads),
    times = 10,
    check = "identical"
)
#> Unit: milliseconds
#>            expr    min     lq  mean median    uq   max neval
#>  cbRcppAlgosSer 357.92 359.18 361.2  360.9 363.6 364.0    10
#>  cbRcppAlgosPar  93.46  96.63 103.1  100.3 101.7 133.7    10
```

## Iterators

We will show one example from each category to demonstrate the
efficiency of the iterators in `RcppAlgos`. The results are similar for
the rest of the cases not shown.

### Combinations

``` r

pkg_arrangements <- function(n, total) {
    a <- icombinations(n, as.integer(n / 2))
    for (i in 1:total) a$getnext()
}

pkg_RcppAlgos <- function(n, total) {
    a <- comboIter(n, as.integer(n / 2))
    for (i in 1:total) a@nextIter()
}

total <- comboCount(18, 9)
total
#> [1] 48620

microbenchmark(cbRcppAlgos    = pkg_RcppAlgos(18, total),
               cbArrangements = pkg_arrangements(18, total),
               times = 15, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median   uq   max neval
#>     cbRcppAlgos  1.00  1.00  1.00   1.00  1.0  1.00    15
#>  cbArrangements 19.52 19.19 19.02  19.12 18.6 19.27    15
```

### Permutations

``` r

pkg_arrangements <- function(n, total) {
    a <- ipermutations(n)
    for (i in 1:total) a$getnext()
}

pkg_RcppAlgos <- function(n, total) {
    a <- permuteIter(n)
    for (i in 1:total) a@nextIter()
}

total <- permuteCount(8)
total
#> [1] 40320

microbenchmark(cbRcppAlgos    = pkg_RcppAlgos(8, total),
               cbArrangements = pkg_arrangements(8, total),
               times = 15, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>     cbRcppAlgos  1.00  1.00  1.00   1.00  1.00  1.00    15
#>  cbArrangements 18.93 18.69 18.59  18.57 18.27 19.37    15
```

### Partitions

``` r

pkg_partitions <- function(n, total) {
    a <- firstpart(n)
    for (i in 1:(total - 1)) a <- nextpart(a)
}

pkg_arrangements <- function(n, total) {
    a <- ipartitions(n)
    for (i in 1:total) a$getnext()
}

pkg_RcppAlgos <- function(n, total) {
    a <- partitionsIter(0:n, repetition = TRUE)
    for (i in 1:total) a@nextIter()
}

total <- partitionsCount(0:40, repetition = TRUE)
total
#> [1] 37338

microbenchmark(cbRcppAlgos    = pkg_RcppAlgos(40, total),
               cbArrangements = pkg_arrangements(40, total),
               cbPartitions   = pkg_partitions(40, total),
               times = 15, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>     cbRcppAlgos  1.00  1.00  1.00   1.00  1.00  1.00    15
#>  cbArrangements 15.29 15.30 14.98  15.37 15.29 13.68    15
#>    cbPartitions 24.63 24.87 24.13  24.98 24.63 21.26    15
```

### Compositions

``` r

pkg_partitions <- function(n, total) {
    a <- firstcomposition(n)
    for (i in 1:(total - 1)) a <- nextcomposition(a, FALSE)
}

pkg_arrangements <- function(n, total) {
    a <- icompositions(n)
    for (i in 1:total) a$getnext()
}

pkg_RcppAlgos <- function(n, total) {
    a <- compositionsIter(0:n, repetition = TRUE)
    for (i in 1:total) a@nextIter()
}

total <- compositionsCount(0:15, repetition = TRUE)
total
#> [1] 16384

microbenchmark(cbRcppAlgos    = pkg_RcppAlgos(15, total),
               cbArrangements = pkg_arrangements(15, total),
               cbPartitions   = pkg_partitions(15, total),
               times = 15, unit = "relative")
#> Unit: relative
#>            expr   min    lq  mean median    uq   max neval
#>     cbRcppAlgos  1.00  1.00  1.00   1.00  1.00  1.00    15
#>  cbArrangements 13.95 13.88 13.70  13.80 13.71 12.21    15
#>    cbPartitions 46.39 45.92 46.04  46.32 46.19 41.90    15
```
