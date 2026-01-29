# Apply Divisor Function to Every Element in a Range

Sieve that generates the number of divisors for every number between
`bound1` and `bound2` (if supplied) or all numbers up to `bound1`. This
is equivalent to applying the divisor function (often written as
\\\sigma(x)\\) to every number in a given range.

## Usage

``` r
numDivisorSieve(bound1, bound2 = NULL, namedVector = FALSE, nThreads = NULL)
```

## Arguments

- bound1:

  Positive integer or numeric value.

- bound2:

  Positive integer or numeric value.

- namedVector:

  Logical flag. If `TRUE`, a named vector is returned. The default is
  `FALSE`.

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

## Details

Simple and efficient sieve that calculates the number of divisors for
every number in a given range. This function is very useful when you
need to calculate the number of divisors for many numbers.

This algorithm benefits greatly from the fast integer division library
'libdivide'. The following is from <https://libdivide.com/>:

- “*libdivide allows you to replace expensive integer divides with
  comparatively cheap multiplication and bitshifts. Compilers usually do
  this, but only when the divisor is known at compile time. libdivide
  allows you to take advantage of it at runtime. The result is that
  integer division can become faster - a lot faster.*”

## Value

Returns a named/unnamed integer vector

## Author

Joseph Wood

## Note

The maximum allowed value is \\2^{53} - 1\\.

## References

- [Divisor function](https://en.wikipedia.org/wiki/Divisor_function)

- [ridiculousfish (author of libdivide)](https://ridiculousfish.com/)

- [github.com/ridiculousfish/libdivide](https://github.com/ridiculousfish/libdivide)

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## Examples

``` r
## Generate some random data
set.seed(8128)
mySamp <- sample(10^6, 5*10^5)

## Generate number of divisors for
## every number less than a million
system.time(mySigmas <- numDivisorSieve(10^6))
#>    user  system elapsed 
#>   0.008   0.000   0.008 

## Now use result in algorithm
for (s in mySamp) {
    sSig <- mySigmas[s]
    ## Continue algorithm
}

## Generating number of divisors for every
## number in a range is no problem
system.time(sigmaRange <- numDivisorSieve(10^13, 10^13 + 10^6))
#>    user  system elapsed 
#>   0.024   0.001   0.024 

## Returning a named vector
numDivisorSieve(10, 20, namedVector = TRUE)
#> 10 11 12 13 14 15 16 17 18 19 20 
#>  4  2  6  2  4  4  5  2  6  2  6 
numDivisorSieve(10, namedVector = TRUE)
#>  1  2  3  4  5  6  7  8  9 10 
#>  1  2  2  3  2  4  2  4  3  4 

## Using nThreads
system.time(numDivisorSieve(1e5, 2e5, nThreads = 2))
#>    user  system elapsed 
#>   0.001   0.000   0.000 
```
