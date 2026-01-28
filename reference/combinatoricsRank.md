# Rank Combinations and Permutations

- Generate the rank (lexicographically) of combinations/permutations.
  These functions are the complement to `comboSample` and
  `permuteSample`. See the examples below.

- GMP support allows for exploration of combinations/permutations of
  vectors with many elements.

## Usage

``` r
comboRank(..., v, repetition = FALSE, freqs = NULL)
permuteRank(..., v, repetition = FALSE, freqs = NULL)
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

  Logical value indicating whether combinations/permutations should be
  with or without repetition. The default is `FALSE`.

- freqs:

  A vector of frequencies used for producing all
  combinations/permutations of a multiset of `v`. Each element of
  `freqs` represents how many times each element of the source vector,
  `v`, is repeated. It is analogous to the `times` argument in
  [`rep`](https://rdrr.io/r/base/rep.html). The default value is `NULL`.

## Details

These algorithms rely on efficiently ranking the \\n^{th}\\
lexicographical combination/permutation.

## Value

A vector of class `integer`, `numeric`, or `bigz` determined by the
total number of combinations/permutations

## References

[Lexicographical
order](https://en.wikipedia.org/wiki/Lexicographical_order)
[ranking/unranking](https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation)

## See also

[`comboSample`](https://jwood000.github.io/RcppAlgos/reference/combinatoricsSample.md),
[`permuteSample`](https://jwood000.github.io/RcppAlgos/reference/combinatoricsSample.md)

## Author

Joseph Wood

## Note

`v` must be supplied.

## Examples

``` r
mySamp = comboSample(30, 8, TRUE, n = 5, seed = 10, namedSample = TRUE)
myRank = comboRank(mySamp, v = 30, repetition = TRUE)
all.equal(as.integer(rownames(mySamp)), myRank)
#> [1] TRUE
```
