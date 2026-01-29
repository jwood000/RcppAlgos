# Generate Complete Factorization for Numbers in a Range

Sieve that generates the complete factorization of all numbers between
`bound1` and `bound2` (if supplied) or all numbers up to `bound1`.

## Usage

``` r
divisorsSieve(bound1, bound2 = NULL, namedList = FALSE, nThreads = NULL)
```

## Arguments

- bound1:

  Positive integer or numeric value.

- bound2:

  Positive integer or numeric value.

- namedList:

  Logical flag. If `TRUE`, a named list is returned. The default is
  `FALSE`.

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

## Details

This function is useful when many complete factorizations are needed.
Instead of generating the complete factorization on the fly, one can
reference the indices/names of the generated list.

This algorithm benefits greatly from the fast integer division library
'libdivide'. The following is from <https://libdivide.com/>:

- “*libdivide allows you to replace expensive integer divides with
  comparatively cheap multiplication and bitshifts. Compilers usually do
  this, but only when the divisor is known at compile time. libdivide
  allows you to take advantage of it at runtime. The result is that
  integer division can become faster - a lot faster.*”

## Value

Returns a named/unnamed list of integer vectors if `max(bound1, bound2)`
\\\< 2^{31}\\, or a list of numeric vectors otherwise.

## Author

Joseph Wood

## Note

The maximum value for either of the bounds is \\2^{53} - 1\\.

## References

- [Divisor](https://en.wikipedia.org/wiki/Divisor)

- [ridiculousfish (author of libdivide)](https://ridiculousfish.com/)

- [github.com/ridiculousfish/libdivide](https://github.com/ridiculousfish/libdivide)

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## See also

[`divisorsRcpp`](https://jwood000.github.io/RcppAlgos/reference/divisorsRcpp.md),
[`primeFactorizeSieve`](https://jwood000.github.io/RcppAlgos/reference/primeFactorizeSieve.md)

## Examples

``` r
## Generate some random data
set.seed(33550336)
mySamp <- sample(10^5, 5*10^4)

## Generate complete factorizations up
## to 10^5 (max element from mySamp)
system.time(allFacs <- divisorsSieve(10^5))
#>    user  system elapsed 
#>   0.024   0.000   0.024 

## Use generated complete factorization for further
## analysis by accessing the index of allFacs
for (s in mySamp) {
    myFac <- allFacs[[s]]
    ## Continue algorithm
}

## Generating complete factorizations over
## a range is efficient as well
system.time(divisorsSieve(10^12, 10^12 + 10^5))
#>    user  system elapsed 
#>   0.060   0.005   0.065 

## Use nThreads for improved efficiency
system.time(divisorsSieve(10^12, 10^12 + 10^5, nThreads = 2))
#>    user  system elapsed 
#>   0.074   0.019   0.059 

## Set 'namedList' to TRUE to return a named list
divisorsSieve(27, 30, namedList = TRUE)
#> $`27`
#> [1]  1  3  9 27
#> 
#> $`28`
#> [1]  1  2  4  7 14 28
#> 
#> $`29`
#> [1]  1 29
#> 
#> $`30`
#> [1]  1  2  3  5  6 10 15 30
#> 

## Using nThreads
system.time(divisorsSieve(1e5, 2e5, nThreads = 2))
#>    user  system elapsed 
#>   0.028   0.001   0.023 
```
