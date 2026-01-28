# Number of Partitions/Compositions

Calculate the number of partitions/compositions of a vector chosen \\m\\
at a time with or without replacement. Additionally, these functions can
calculate the number of partitions of multisets.

## Usage

``` r
partitionsCount(v, m = NULL, ...)
compositionsCount(v, m = NULL, ...)

# Default S3 method
partitionsCount(v, m = NULL, repetition = FALSE,
                freqs = NULL, target = NULL, ...)
# Default S3 method
compositionsCount(v, m = NULL, repetition = FALSE,
                  freqs = NULL, target = NULL, weak = FALSE, ...)

# S3 method for class 'table'
partitionsCount(v, m = NULL, target = NULL, ...)
# S3 method for class 'table'
compositionsCount(v, m = NULL, target = NULL, weak = FALSE, ...)
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

## Value

A numerical value representing the total number of
partitions/compositions.

## Note

When the number of results exceeds \\2^{53} - 1\\, a number of class
`bigz` is returned.

## See also

[`partitionsGeneral`](https://jwood000.github.io/RcppAlgos/reference/partitionsGeneral.md),
[`compositionsGeneral`](https://jwood000.github.io/RcppAlgos/reference/partitionsGeneral.md)

## Examples

``` r
## Same interface as partitionsGeneral
partitionsCount(25, 5)
#> [1] 30
compositionsCount(25, 5, TRUE)
#> [1] 10626
partitionsCount(15, 7, TRUE)
#> [1] 21
partitionsCount(25, 5, freqs = rep(2, 25))
#> [1] 148

## Return object of class 'bigz'
partitionsCount(2500, 15, TRUE)
#> Big Integer ('bigz') :
#> [1] 4188807783906550542066742
compositionsCount(2500, 15, TRUE)
#> Big Integer ('bigz') :
#> [1] 4097094997743241086270928785678827376
```
