# Apply Euler's Phi Function to Every Element in a Range

Sieve that generates the number of coprime elements for every number
between `bound1` and `bound2` (if supplied) or all numbers up to
`bound1`. This is equivalent to applying Euler's phi function (often
written as \\\phi(x)\\) to every number in a given range.

## Usage

``` r
eulerPhiSieve(bound1, bound2 = NULL, namedVector = FALSE, nThreads = NULL)
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

For the simple case (i.e. when `bound2 = NULL`), this algorithm first
generates all primes up to \\n\\ via the sieve of Eratosthenes. We use
these primes to sieve over the sequence `1:n`, dividing each value by
\\p\\, creating a temporary value that will be subtracted from the
original value at each index (i.e. equivalent to multiply each index by
\\(1 - 1/p)\\ but more efficient as we don't have to deal with floating
point numbers). The case when `is.null(bound2) = FALSE` is more
complicated but the basic ideas still hold.

This function is very useful when you need to calculate Euler's phi
function for many numbers in a range as performing this calculation on
the fly can be computationally expensive.

This algorithm benefits greatly from the fast integer division library
'libdivide'. The following is from <https://libdivide.com/>:

- “*libdivide allows you to replace expensive integer divides with
  comparatively cheap multiplication and bitshifts. Compilers usually do
  this, but only when the divisor is known at compile time. libdivide
  allows you to take advantage of it at runtime. The result is that
  integer division can become faster - a lot faster.*”

## Value

Returns a named/unnamed integer vector if `max(bound1, bound2)` \\\<
2^{31}\\, or a numeric vector otherwise.

## Author

Joseph Wood

## Note

The maximum allowed value is \\2^{53} - 1\\.

## References

- [Euler's totient
  function](https://en.wikipedia.org/wiki/Euler%27s_totient_function)

- [ridiculousfish (author of libdivide)](https://ridiculousfish.com/)

- [github.com/ridiculousfish/libdivide](https://github.com/ridiculousfish/libdivide)

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## Examples

``` r
## Generate some random data
set.seed(496)
mySamp <- sample(10^6, 5*10^5)

## Generate number of coprime elements for many numbers
system.time(myPhis <- eulerPhiSieve(10^6))
#>    user  system elapsed 
#>   0.009   0.001   0.009 

## Now use result in algorithm
for (s in mySamp) {
    sPhi <- myPhis[s]
    ## Continue algorithm
}

## See https://projecteuler.net
system.time(which.max((1:10^6)/eulerPhiSieve(10^6)))
#>    user  system elapsed 
#>   0.011   0.001   0.012 

## Generating number of coprime elements
## for every number in a range is no problem
system.time(myPhiRange <- eulerPhiSieve(10^13, 10^13 + 10^6))
#>    user  system elapsed 
#>   0.027   0.002   0.029 

## Returning a named vector
eulerPhiSieve(10, 20, namedVector = TRUE)
#> 10 11 12 13 14 15 16 17 18 19 20 
#>  4 10  4 12  6  8  8 16  6 18  8 
eulerPhiSieve(10, namedVector = TRUE)
#>  1  2  3  4  5  6  7  8  9 10 
#>  1  1  2  2  4  2  6  4  6  4 

## Using nThreads
system.time(eulerPhiSieve(1e5, 2e5, nThreads = 2))
#>    user  system elapsed 
#>   0.001   0.000   0.001 
```
