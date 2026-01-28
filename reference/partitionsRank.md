# Rank Partitions/Compositions

- Generate the rank (lexicographically) of partitions/compositions.
  These functions are the complement to `partitions/compositionsSample`.
  See the examples below.

- GMP support allows for exploration of partitions/compositions of
  vectors with many elements.

## Usage

``` r
partitionsRank(..., v, repetition = FALSE, freqs = NULL, target = NULL)

compositionsRank(..., v, repetition = FALSE, freqs = NULL,
                 target = NULL, weak = FALSE)
```

## Arguments

- ...:

  vectors or matrices to be ranked.

- v:

  Source vector. If `v` is a positive integer, it will be converted to
  the sequence `1:v`. If `v` is a negative integer, it will be converted
  to the sequence `v:-1`. All atomic types are supported (See
  [`is.atomic`](https://rdrr.io/r/base/is.recursive.html)).

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

## Details

These algorithms rely on efficiently ranking the \\n^{th}\\
lexicographical partition.

## Value

A vector of class `integer`, `numeric`, or `bigz` determined by the
total number of partitions/compositions

## References

[Lexicographical
order](https://en.wikipedia.org/wiki/Lexicographical_order)
[ranking/unranking](https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation)

## See also

[`partitionsSample`](https://jwood000.github.io/RcppAlgos/reference/partitionsSample.md),
[`compositionsSample`](https://jwood000.github.io/RcppAlgos/reference/partitionsSample.md)

## Author

Joseph Wood

## Note

`v` must be supplied.

## Examples

``` r
mySamp = partitionsSample(30, 8, TRUE, n = 5, seed = 10, namedSample = TRUE)
myRank = partitionsRank(mySamp, v = 30, repetition = TRUE)
all.equal(as.integer(rownames(mySamp)), myRank)
#> [1] TRUE
```
