# Vectorized Factorization (Complete)

Function for generating the complete factorization for a vector of
numbers.

## Usage

``` r
divisorsRcpp(v, namedList = FALSE, nThreads = NULL)
```

## Arguments

- v:

  Vector of integers or numeric values. Non-integral values will be
  coerced to whole numbers.

- namedList:

  Logical flag. If `TRUE` and the `length(v) > 1`, a named list is
  returned. The default is `FALSE`.

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

## Details

Efficient algorithm that builds on
[`primeFactorize`](https://jwood000.github.io/RcppAlgos/reference/primeFactorize.md)
to generate the complete factorization of many numbers.

## Value

- Returns an unnamed vector if `length(v) == 1` regardless of the value
  of `namedList`. If \\v \< 2^{31}\\, the class of the returned vector
  will be integer, otherwise the class will be numeric.

- If `length(v) > 1`, a named/unnamed list of vectors will be returned.
  If `max(bound1, bound2)` \\\< 2^{31}\\, the class of each vector will
  be integer, otherwise the class will be numeric.

## References

- [Divisor](https://en.wikipedia.org/wiki/Divisor)

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## Author

Joseph Wood

## Note

The maximum value for each element in \\v\\ is \\2^{53} - 1\\.

## See also

[`primeFactorize`](https://jwood000.github.io/RcppAlgos/reference/primeFactorize.md)

## Examples

``` r
## Get the complete factorization of a single number
divisorsRcpp(10^8)
#>  [1]         1         2         4         5         8        10        16
#>  [8]        20        25        32        40        50        64        80
#> [15]       100       125       128       160       200       250       256
#> [22]       320       400       500       625       640       800      1000
#> [29]      1250      1280      1600      2000      2500      3125      3200
#> [36]      4000      5000      6250      6400      8000     10000     12500
#> [43]     15625     16000     20000     25000     31250     32000     40000
#> [50]     50000     62500     78125     80000    100000    125000    156250
#> [57]    160000    200000    250000    312500    390625    400000    500000
#> [64]    625000    781250    800000   1000000   1250000   1562500   2000000
#> [71]   2500000   3125000   4000000   5000000   6250000  10000000  12500000
#> [78]  20000000  25000000  50000000 100000000

## Or get the complete factorization of many numbers
set.seed(29)
myVec <- sample(-1000000:1000000, 1000)
system.time(myFacs <- divisorsRcpp(myVec))
#>    user  system elapsed 
#>   0.001   0.000   0.001 

## Return named list
myFacsWithNames <- divisorsRcpp(myVec, namedList = TRUE)

## Using nThreads
system.time(divisorsRcpp(myVec, nThreads = 2))
#>    user  system elapsed 
#>   0.002   0.000   0.001 
```
