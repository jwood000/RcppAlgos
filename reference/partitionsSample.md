# Sample Partitions/Compositions

- Generate a specific (lexicographically) or random sample of
  partitions/compositions of a number.

- Produce results in parallel using the `Parallel` or `nThreads`
  arguments.

- GMP support allows for exploration of cases where the number of
  partitions/compositions is large.

## Usage

``` r
partitionsSample(v, m = NULL, ...)
compositionsSample(v, m = NULL, ...)

# Default S3 method
partitionsSample(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL,
    n = NULL, sampleVec = NULL, seed = NULL,
    nThreads = NULL, namedSample = FALSE, ...
)
# Default S3 method
compositionsSample(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL,
    weak = FALSE, n = NULL, sampleVec = NULL, seed = NULL,
    nThreads = NULL, namedSample = FALSE, ...
)

# S3 method for class 'table'
partitionsSample(
    v, m = NULL, target = NULL, n = NULL,
    sampleVec = NULL, seed = NULL, nThreads = NULL, namedSample = FALSE, ...
)
# S3 method for class 'table'
compositionsSample(
    v, m = NULL, target = NULL, weak = FALSE, n = NULL,
    sampleVec = NULL, seed = NULL, nThreads = NULL, namedSample = FALSE, ...
)
```

## Arguments

- v:

  Source vector. If `v` is a positive integer, it will be converted to
  the sequence `1:v`. If `v` is a negative integer, it will be converted
  to the sequence `v:-1`. Only integer and numeric vectors are accepted.

- m:

  Width of the partition. If `m = NULL`, the length will be determined
  by the partitioning case (*e.g.* When we are generating distinct
  partitions of \\n\\, the width will be equal to the smallest \\m\\
  such that `sum(1:m) >= n`).

- ...:

  Further arguments passed to methods.

- repetition:

  Logical value indicating whether partitions/compositions should be
  with or without repetition. The default is `FALSE`.

- freqs:

  A vector of frequencies used for producing all partitions of a
  multiset of `v`. Each element of `freqs` represents how many times
  each element of the source vector, `v`, is repeated. It is analogous
  to the `times` argument in [`rep`](https://rdrr.io/r/base/rep.html).
  The default value is `NULL`.

- target:

  Number to be partitioned. If `NULL`, `max(v)` will be used.

- weak:

  (Compositions only) Logical flag indicating whether to allow terms of
  the sequence to be zero.

- n:

  Number of partitions/compositions to return. The default is `NULL`.

- sampleVec:

  A vector of numbers representing the lexicographical
  partitions/compositions to return. Accepts vectors of class `bigz` as
  well as vectors of characters

- seed:

  Random seed initialization. The default is `NULL`. N.B. If the gmp
  library is needed, this parameter must be set in order to have
  reproducible results (*E.g*
  [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect in
  these cases).

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

- namedSample:

  Logical flag. If `TRUE`, `rownames` corresponding to the
  lexicographical partition, will be added to the returned matrix. The
  default is `FALSE`.

## Details

These algorithms rely on efficiently generating the \\n^{th}\\
lexicographical partition. This is the process of
[unranking](https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation).

## Value

A matrix is returned with each row containing a vector of length \\m\\.

## References

- [Lexicographical
  order](https://en.wikipedia.org/wiki/Lexicographical_order)

- [Partition (Number
  Theory)](https://en.wikipedia.org/wiki/Partition_(number_theory))

## Author

Joseph Wood

## Note

- `partitionsSample` is not available for the following cases:

  - With standard multisets. If zero is the only element with a
    non-trivial multiplicity, sampling is allowed (*e.g.*
    `partitionsSample(0:100, freqs = c(100, rep(1, 100)), n = 2)`)

  - If the source vector is not isomorphic to `1:length(v)` (*e.g.*
    `v = c(1, 4, 6, 7, 8)`).

- `n` and `sampleVec` cannot both be `NULL`.

## Examples

``` r
partitionsSample(100, 10, n = 5)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    2    3    4    6    8    9   13   14   17    24
#> [2,]    2    3    4    6    8    9   11   12   20    25
#> [3,]    1    2    6    7    8   10   12   15   18    21
#> [4,]    1    3    6    7    9   11   12   16   17    18
#> [5,]    1    2    5    6    7    8   10   15   18    28
partitionsSample(100, 10, seed = 42, n = 5, target = 200)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    1    6    9   11   12   14   20   25   48    54
#> [2,]    1    3    6   11   12   13   17   26   49    62
#> [3,]    3    5    6    7    8   21   24   36   41    49
#> [4,]    3    6    8    9   22   27   28   29   33    35
#> [5,]    1    2    3    8   13   29   33   34   36    41

## retrieve specific results (lexicographically)
partitionsCount(100, 10, TRUE, target = 500)
#> [1] 175591757896
## [1] 175591757896
partitionsSample(100, 10, TRUE, target = 500,
                 sampleVec = c(1, 1000, 175591757896))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    1    1    1    1    1   95  100  100  100   100
#> [2,]    1    1    1    1   16   90   94   96  100   100
#> [3,]   50   50   50   50   50   50   50   50   50    50
```
