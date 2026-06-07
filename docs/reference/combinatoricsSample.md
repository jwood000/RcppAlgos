# Sample Combinations and Permutations

- Generate a specific (lexicographically) or random sample of
  combinations/permutations.

- Produce results in parallel using the `Parallel` or `nThreads`
  arguments.

- GMP support allows for exploration of combinations/permutations of
  vectors with many elements.

## Usage

``` r
comboSample(v, m = NULL, ...)
permuteSample(v, m = NULL, ...)

# S3 method for class 'numeric'
comboSample(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
            sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
            nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...)

# S3 method for class 'numeric'
permuteSample(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
              sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
              nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...)

# S3 method for class 'factor'
comboSample(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
    nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
)
# S3 method for class 'factor'
permuteSample(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
    nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
)

# Default S3 method
comboSample(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL, sampleVec = NULL,
    seed = NULL, FUN = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
)
# Default S3 method
permuteSample(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL, sampleVec = NULL,
    seed = NULL, FUN = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
)

# S3 method for class 'table'
comboSample(
    v, m = NULL, n = NULL, sampleVec = NULL, seed = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
)
# S3 method for class 'table'
permuteSample(
    v, m = NULL, n = NULL, sampleVec = NULL, seed = NULL, FUN = NULL,
    Parallel = FALSE, nThreads = NULL, namedSample = FALSE, FUN.VALUE = NULL, ...
)

# S3 method for class 'list'
comboSample(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, namedSample = FALSE, ...
)
# S3 method for class 'list'
permuteSample(
    v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, namedSample = FALSE, ...
)
```

## Arguments

- v:

  Source vector. If `v` is a positive integer, it will be converted to
  the sequence `1:v`. If `v` is a negative integer, it will be converted
  to the sequence `v:-1`. All atomic types are supported (See
  [`is.atomic`](https://rdrr.io/r/base/is.recursive.html)).

- m:

  Number of elements to choose. If `repetition = TRUE` or `freqs` is
  utilized, `m` can exceed the length of `v`. If `m = NULL`, the length
  will default to `length(v)` or `sum(freqs)`.

- ...:

  Further arguments passed to methods.

- repetition:

  Logical value indicating whether combinations/permutations should be
  with or without repetition. The default is `FALSE`.

- freqs:

  A vector of frequencies used for producing all
  combinations/permutations of a multiset of `v`. Each element of
  `freqs` represents how many times each element of the source vector,
  `v`, is repeated. It is analogous to the `times` argument in
  [`rep`](https://rdrr.io/r/base/rep.html). The default value is `NULL`.

- n:

  Number of combinations/permutations to return. The default is `NULL`.

- sampleVec:

  A vector of indices representing the lexicographical
  combination/permutations to return. Accepts whole numbers as well as
  vectors of class `bigz` as well as vectors of characters

- seed:

  Random seed initialization. The default is `NULL`. N.B. If the gmp
  library is needed, this parameter must be set in order to have
  reproducible results (*E.g*
  [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect in
  these cases).

- FUN:

  Function to be applied to each combination/permutation. The default is
  `NULL`.

- Parallel:

  Logical value indicating whether combinations/permutations should be
  generated in parallel. The default is `FALSE`. If `TRUE` and
  `nThreads = NULL`, the number of threads used is equal to the minimum
  of one minus the number of threads available on your system and the
  number of results requested (*e.g.* if user has 16 threads and only
  needs 5 results, 5 threads will be used (*i.e.*
  `min(16 - 1, 5) = 5`)). If `nThreads` is not `NULL`, it will be given
  preference (*e.g.* if user has 8 threads with `Parallel = TRUE` and
  `nThreads = 4`, only 4 threads will be spawned). If your system is
  single-threaded, the arguments `Parallel` and `nThreads` are ignored.

- nThreads:

  Specific number of threads to be used. The default is `NULL`. See
  `Parallel`.

- namedSample:

  Logical flag. If `TRUE`, `rownames` corresponding to the
  lexicographical combination/permutation, will be added to the returned
  matrix. The default is `FALSE`.

- FUN.VALUE:

  A template for the return value from `FUN`. See 'Details' of
  [`vapply`](https://rdrr.io/r/base/lapply.html) for more information.

## Details

These algorithms rely on efficiently generating the \\n^{th}\\
lexicographical combination/permutation. This is the process of
[unranking](https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation).

## Value

- In general, a matrix with \\m\\ or \\m + 1\\ columns, depending on the
  value of `keepResults`

- If `FUN` is utilized and `FUN.VALUE = NULL`, a list is returned

- When both `FUN` and `FUN.VALUE` are not `NULL`, the return is modeled
  after the return of `vapply`. See the 'Value' section of
  [`vapply`](https://rdrr.io/r/base/lapply.html).

## References

[Lexicographical
order](https://en.wikipedia.org/wiki/Lexicographical_order)

## See also

[`comboRank`](https://jwood000.github.io/RcppAlgos/reference/combinatoricsRank.md),
[`permuteRank`](https://jwood000.github.io/RcppAlgos/reference/combinatoricsRank.md)

## Author

Joseph Wood

## Note

- `Parallel` and `nThreads` will be ignored in the following cases:

  - If the class of the vector passed is `character` (N.B.
    `Rcpp::CharacterMatrix` is not thread safe). Alternatively, you can
    generate an indexing matrix in parallel.

  - If `FUN` is utilized.

- `n` and `sampleVec` cannot both be `NULL`.

- Factor vectors are accepted. Class and level attributes are preserved
  except when `FUN` is used.

## Examples

``` r
## generate 10 random combinations
comboSample(30, 8, TRUE, n = 5, seed = 10)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    7   11   11   17   20   22   23   27
#> [2,]    4    5    6   15   20   22   22   28
#> [3,]    6   11   11   17   17   21   24   24
#> [4,]   15   17   18   25   27   27   29   30
#> [5,]    5   12   14   18   21   22   25   27

## Using sampleVec to generate specific permutations
fqs   = c(1,2,2,1,2,2,1,2,1,2,2,1,2,1,1)
s_idx = c(1, 10^2, 10^5, 10^8, 10^11)

permuteSample(15, 10, freqs = fqs, sampleVec = s_idx)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    1    2    2    3    3    4    5    5    6     6
#> [2,]    1    2    2    3    3    4    5    6    5    10
#> [3,]    1    2    2    3    3   11    8    4   15     6
#> [4,]    1    2    5    3   10    7    6   11    9    13
#> [5,]   13   15    8   11    5    5    6    1   11    12

## Same example using 'table' method
permuteSample(table(rep(1:15, times = fqs)), 10, sampleVec = s_idx)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    1    2    2    3    3    4    5    5    6     6
#> [2,]    1    2    2    3    3    4    5    6    5    10
#> [3,]    1    2    2    3    3   11    8    4   15     6
#> [4,]    1    2    5    3   10    7    6   11    9    13
#> [5,]   13   15    8   11    5    5    6    1   11    12

## Generate each result one by one...
## Same, but not as efficient as generating iteratively
all.equal(comboSample(10, 5, sampleVec = 1:comboCount(10, 5)),
          comboGeneral(10, 5))
#> [1] TRUE

## Examples with enormous number of total permutations
num = permuteCount(10000, 20)
gmp::log2.bigz(num)
#> [1] 265.7268

first  = gmp::urand.bigz(n = 1, size = 265, seed = 123)
#> Seed default initialisation
#> Seed initialisation
mySamp = do.call(c, lapply(0:10, function(x) gmp::add.bigz(first, x)))

class(mySamp)
#> [1] "bigz"

## using permuteSample
pSamp = permuteSample(10000, 20, sampleVec = mySamp)

## using permuteGeneral
pGeneral = permuteGeneral(10000, 20,
                          lower = first,
                          upper = gmp::add.bigz(first, 10))

identical(pSamp, pGeneral)
#> [1] TRUE

## Using nThreads
permPar = permuteSample(10000, 50, n = 8, seed = 10, nThreads = 2)

## Using FUN
permuteSample(10000, 50, n = 4, seed = 10, FUN = sd)
#> [[1]]
#> [1] 2694.073
#> 
#> [[2]]
#> [1] 2923.153
#> 
#> [[3]]
#> [1] 2971.491
#> 
#> [[4]]
#> [1] 2711.457
#> 

if (FALSE) { # \dontrun{
## Using Parallel
permuteSample(10000, 50, n = 80, seed = 10, Parallel = TRUE)
} # }
```
