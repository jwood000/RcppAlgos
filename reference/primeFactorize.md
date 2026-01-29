# Vectorized Prime Factorization

Implementation of Pollard's rho algorithm for generating the prime
factorization. The algorithm is based on the "factorize.c" source file
from the gmp library found here <https://gmplib.org>.

## Usage

``` r
primeFactorize(v, namedList = FALSE, nThreads = NULL)
```

## Arguments

- v:

  Vector of integers or numeric values. Non-integral values will be
  cured to whole numbers.

- namedList:

  Logical flag. If `TRUE` and the `length(v) > 1`, a named list is
  returned. The default is `FALSE`.

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

## Details

As noted in the Description section above, this algorithm is based on
the "factorize.c" source code from the gmp library. Much of the code in
RcppAlgos::primeFactorize is a straightforward translation from multiple
precision C data types to standard C++ data types. A crucial part of the
algorithm's efficiency is based on quickly determining
[primality](https://en.wikipedia.org/wiki/Primality_test), which is
easily computed with gmp. However, with standard C++, this is quite
challenging. Much of the research for RcppAlgos::primeFactorize was
focused on developing an algorithm that could accurately and efficiently
compute primality.

For more details, see the documentation for
[`isPrimeRcpp`](https://jwood000.github.io/RcppAlgos/reference/isPrimeRcpp.md).

## Value

- Returns an unnamed vector if `length(v) == 1` regardless of the value
  of `namedList`. If \\v \< 2^{31}\\, the class of the returned vector
  will be integer, otherwise the class will be numeric.

- If `length(v) > 1`, a named/unnamed list of vectors will be returned.
  If `max(bound1, bound2)` \\\< 2^{31}\\, the class of each vector will
  be integer, otherwise the class will be numeric.

## References

- [Pollard's rho
  algorithm](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm)

- [Miller-Rabin primality
  test](https://en.wikipedia.org/wiki/Miller-Rabin_primality_test)

- [Accurate Modular Arithmetic with Double
  Precision](https://codereview.stackexchange.com/questions/186751/accurate-modular-arithmetic-with-double-precision)

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## Author

Joseph Wood

## Note

The maximum value for each element in \\v\\ is \\2^{53} - 1\\.

## See also

[`primeFactorizeSieve`](https://jwood000.github.io/RcppAlgos/reference/primeFactorizeSieve.md),
[`factorize`](https://rdrr.io/pkg/gmp/man/factor.html)

## Examples

``` r
## Get the prime factorization of a single number
primeFactorize(10^8)
#>  [1] 2 2 2 2 2 2 2 2 5 5 5 5 5 5 5 5

## Or get the prime factorization of many numbers
set.seed(29)
myVec <- sample(-1000000:1000000, 1000)
system.time(pFacs <- primeFactorize(myVec))
#>    user  system elapsed 
#>   0.001   0.000   0.001 

## Return named list
pFacsWithNames <- primeFactorize(myVec, namedList = TRUE)

## Using nThreads
system.time(primeFactorize(myVec, nThreads = 2))
#>    user  system elapsed 
#>       0       0       0 
```
