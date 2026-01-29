# Generate Prime Factorization for Numbers in a Range

Generates the prime factorization of all numbers between `bound1` and
`bound2` (if supplied) or all numbers up to `bound1`.

## Usage

``` r
primeFactorizeSieve(bound1, bound2 = NULL, namedList = FALSE, nThreads = NULL)
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

This function is useful when many prime factorizations are needed.
Instead of generating the prime factorization on the fly, one can
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

## Note

The maximum value for either of the bounds is \\2^{53} - 1\\.

## References

- [Prime Factor](https://en.wikipedia.org/wiki/Prime_factor)

- [ridiculousfish (author of libdivide)](https://ridiculousfish.com/)

- [github.com/ridiculousfish/libdivide](https://github.com/ridiculousfish/libdivide)

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## Author

Joseph Wood

## See also

[`primeFactorize`](https://jwood000.github.io/RcppAlgos/reference/primeFactorize.md),
[`divisorsSieve`](https://jwood000.github.io/RcppAlgos/reference/divisorsSieve.md),
[`factorize`](https://rdrr.io/pkg/gmp/man/factor.html)

## Examples

``` r
## Generate some random data
set.seed(28)
mySamp <- sample(10^5, 5*10^4)

## Generate prime factorizations up
## to 10^5 (max element from mySamp)
system.time(allPFacs <- primeFactorizeSieve(10^5))
#>    user  system elapsed 
#>    0.02    0.00    0.02 

## Use generated prime factorization for further
## analysis by accessing the index of allPFacs
for (s in mySamp) {
    pFac <- allPFacs[[s]]
    ## Continue algorithm
}

## Generating prime factorizations over
## a range is efficient as well
system.time(primeFactorizeSieve(10^12, 10^12 + 10^5))
#>    user  system elapsed 
#>   0.033   0.000   0.033 

## Set 'namedList' to TRUE to return a named list
primeFactorizeSieve(27, 30, namedList = TRUE)
#> $`27`
#> [1] 3 3 3
#> 
#> $`28`
#> [1] 2 2 7
#> 
#> $`29`
#> [1] 29
#> 
#> $`30`
#> [1] 2 3 5
#> 

## Using nThreads
system.time(primeFactorizeSieve(1e4, 5e4, nThreads = 2))
#>    user  system elapsed 
#>   0.007   0.000   0.006 
```
