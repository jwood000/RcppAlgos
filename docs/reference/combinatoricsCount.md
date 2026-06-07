# Number of combinations/permutations

Calculate the number of combinations/permutations of a vector chosen
\\m\\ at a time with or without replacement. Additionally, these
functions can calculate the number of combinations/permutations of
multisets.

## Usage

``` r
comboCount(v, m = NULL, ...)
permuteCount(v, m = NULL, ...)

# Default S3 method
comboCount(v, m = NULL, repetition = FALSE, freqs = NULL, ...)
# Default S3 method
permuteCount(v, m = NULL, repetition = FALSE, freqs = NULL, ...)

# S3 method for class 'table'
comboCount(v, m = NULL, ...)
# S3 method for class 'table'
permuteCount(v, m = NULL, ...)
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

## Value

A numerical value representing the total number of
combinations/permutations.

## Note

When the number of results exceeds \\2^{53} - 1\\, a number of class
`bigz` is returned.

## See also

[`comboGeneral`](https://jwood000.github.io/RcppAlgos/reference/combinatoricsGeneral.md),
[`permuteGeneral`](https://jwood000.github.io/RcppAlgos/reference/combinatoricsGeneral.md)

## Examples

``` r
## Same interface as the respective "general" functions:
## i.e. comboGeneral & permuteGeneral

permuteCount(-5)
#> [1] 120
permuteCount(5)
#> [1] 120
comboCount(25, 12)
#> [1] 5200300
permuteCount(15, 7, TRUE)
#> [1] 170859375
comboCount(25, 12, freqs = rep(2, 25))
#> [1] 458917850

## Return object of class 'bigz'
comboCount(250, 15, freqs = rep(2, 250))
#> Big Integer ('bigz') :
#> [1] 1035444613157684247678300
```
