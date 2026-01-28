# Number of Partitions of a Vector into Groups

Calculate the number of partitions of a vector into groups. See the
related integer sequences A025035-A025042 at [OEIS](https://oeis.org)
(E.g. [A025036](https://oeis.org/A025036) for Number of partitions of \\
1, 2, ..., 4n \\ into sets of size 4.)

## Usage

``` r
comboGroupsCount(v, numGroups = NULL, grpSizes = NULL)
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

## Value

A numerical value representing the total number of partitions of groups.

## References

[OEIS Integer Sequence A025036](https://oeis.org/A025036)

## Author

Joseph Wood

## Note

When the number of results exceeds \\2^{53} - 1\\, a number of class
`bigz` is returned.

## Examples

``` r
comboGroupsCount(16, 4)
#> [1] 2627625
comboGroupsCount(16, grpSizes = c(1:4, 6))
#> [1] 100900800
comboGroupsCount(28, grpSizes = rep(2:5, each = 2))
#> Big Integer ('bigz') :
#> [1] 15954139019358540000
```
