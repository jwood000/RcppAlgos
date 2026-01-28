# Sample Partitions of a Vector into Groups

- Generate a specific (lexicographically) or random sample of partitions
  of groups.

- Produce results in parallel using the `Parallel` or `nThreads`
  arguments.

- GMP support allows for exploration where the number of results is
  large.

## Usage

``` r
comboGroupsSample(v, numGroups = NULL, grpSizes = NULL, retType = "matrix",
                  n = NULL, sampleVec = NULL, seed = NULL, Parallel = FALSE,
                  nThreads = NULL, namedSample = FALSE)
```

## Arguments

- v:

  Source vector. If `v` is a positive integer, it will be converted to
  the sequence `1:v`. If `v` is a negative integer, it will be converted
  to the sequence `v:-1`. All atomic types are supported (See
  [`is.atomic`](https://rdrr.io/r/base/is.recursive.html)).

- numGroups:

  An Integer. The number of groups that the vector will be partitioned
  into. The default is `NULL`. If provided and `grpSize` is `NULL`, it
  must divide the length of v (if v is a vector) or v (if v is a
  scalar).

- grpSizes:

  A vector of whole numbers representing the size of each group. The
  default is `NULL`. If provided, the sum of the elements must total the
  length of v (if v is a vector) or v (if v is a scalar).

- retType:

  A string, "3Darray" or "matrix", that determines the shape of the
  output. The default is "matrix". Note, "3Darray" can only be used when
  the size of each group is uniform. When the size of each group varies,
  the return output will always be a matrix.

- n:

  Number of results to return. The default is `NULL`.

- sampleVec:

  A vector of numbers representing the lexicographical partition of
  groups to return. Accepts vectors of class `bigz` as well as vectors
  of characters

- seed:

  Random seed initialization. The default is `NULL`. N.B. If the gmp
  library is needed, this parameter must be set in order to have
  reproducible results (*E.g*
  [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect in
  these cases).

- Parallel:

  Logical value indicating whether results should be generated in
  parallel. The default is `FALSE`. If `TRUE` and `nThreads = NULL`, the
  number of threads used is equal to the minimum of one minus the number
  of threads available on your system and the number of results
  requested (*e.g.* if user has 16 threads and only needs 5 results, 5
  threads will be used (*i.e.* `min(16 - 1, 5) = 5`)). If `nThreads` is
  not `NULL`, it will be given preference (*e.g.* if user has 8 threads
  with `Parallel = TRUE` and `nThreads = 4`, only 4 threads will be
  spawned). If your system is single-threaded, the arguments `Parallel`
  and `nThreads` are ignored.

- nThreads:

  Specific number of threads to be used. The default is `NULL`. See
  `Parallel`.

- namedSample:

  Logical flag. If `TRUE`, `rownames` corresponding to the
  lexicographical result, will be added to the returned matrix. The
  default is `FALSE`.

## Details

These algorithms rely on efficiently generating the \\n^{th}\\
lexicographical result.

## Value

By default, a matrix is returned with column names corresponding to the
associated group. If `retType = "3Darray"`, a 3D array is returned.

## References

[Lexicographical
order](https://en.wikipedia.org/wiki/Lexicographical_order)

## Author

Joseph Wood

## Examples

``` r
## generate 10 random partitions of groups of equal size
comboGroupsSample(10, 2, n = 10, seed = 123)
#>       Grp1 Grp1 Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp2 Grp2
#>  [1,]    1    2    4    7    8    3    5    6    9   10
#>  [2,]    1    3    5    8    9    2    4    6    7   10
#>  [3,]    1    2    6    8   10    3    4    5    7    9
#>  [4,]    1    2    3    6    9    4    5    7    8   10
#>  [5,]    1    3    4    7    9    2    5    6    8   10
#>  [6,]    1    2    5    7    9    3    4    6    8   10
#>  [7,]    1    2    6    8    9    3    4    5    7   10
#>  [8,]    1    5    7    8    9    2    3    4    6   10
#>  [9,]    1    2    5    7   10    3    4    6    8    9
#> [10,]    1    4    5    9   10    2    3    6    7    8

## generate 10 random partitions of groups of varying sizes
comboGroupsSample(10, grpSizes = 1:4, n = 10, seed = 123)
#>       Grp1 Grp2 Grp2 Grp3 Grp3 Grp3 Grp4 Grp4 Grp4 Grp4
#>  [1,]    2    8   10    1    6    7    3    4    5    9
#>  [2,]    2    9   10    4    5    6    1    3    7    8
#>  [3,]    9    2    4    3    7   10    1    5    6    8
#>  [4,]    7    8    9    1    2    5    3    4    6   10
#>  [5,]   10    6    9    2    5    7    1    3    4    8
#>  [6,]    3    2    9    1    6    8    4    5    7   10
#>  [7,]    2    4    6    3    7   10    1    5    8    9
#>  [8,]    8    2   10    3    6    9    1    4    5    7
#>  [9,]    3    5    9    1    6    8    2    4    7   10
#> [10,]   10    2    3    4    5    8    1    6    7    9

## using sampleVec to generate specific results
comboGroupsSample(15, 5, sampleVec = c(1, 100, 1e3, 1e6))
#>      Grp1 Grp1 Grp1 Grp2 Grp2 Grp2 Grp3 Grp3 Grp3 Grp4 Grp4 Grp4 Grp5 Grp5 Grp5
#> [1,]    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
#> [2,]    1    2    3    4    5    6    7    9   12    8   14   15   10   11   13
#> [3,]    1    2    3    4    5    9    6   10   13    7   14   15    8   11   12
#> [4,]    1    8   10    2   12   15    3    5   13    4   11   14    6    7    9

all.equal(comboGroupsSample(10, 5,
            sampleVec = 1:comboGroupsCount(10, 5)),
         comboGroups(10, 5))
#> [1] TRUE

## Examples with enormous number of total results
num = comboGroupsCount(100, 20)
gmp::log2.bigz(num)
#> [1] 325.5498
## [1] 325.5498

first = gmp::urand.bigz(n = 1, size = 325, seed = 123)
#> Seed initialisation
mySamp = do.call(c, lapply(0:10, function(x) gmp::add.bigz(first, x)))

class(mySamp)
#> [1] "bigz"
## [1] "bigz"

## using the sampling function
cbgSamp = comboGroupsSample(100, 20, sampleVec = mySamp)

## using the standard function
cbgGeneral = comboGroups(100, 20,
                         lower = first,
                         upper = gmp::add.bigz(first, 10))

identical(cbgSamp, cbgGeneral)
#> [1] TRUE
## [1] TRUE

if (FALSE) { # \dontrun{
## Using Parallel
system.time(comboGroupsSample(1000, 20, n = 80, seed = 10, Parallel = TRUE))
} # }
```
