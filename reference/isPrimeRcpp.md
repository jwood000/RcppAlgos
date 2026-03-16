# Vectorized Primality Test

Implementation of the [Miller-Rabin primality
test](https://en.wikipedia.org/wiki/Miller-Rabin_primality_test). Based
on the "mp_prime_p" function from the "factorize.c" source file found in
the gmp library: <https://gmplib.org>.

## Usage

``` r
isPrimeRcpp(v, namedVector = FALSE, nThreads = NULL)
```

## Arguments

- v:

  Vector of integers or numeric values.

- namedVector:

  Logical flag. If `TRUE`, a named vector is returned. The default is
  `FALSE`.

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

## Note

The maximum value for each element in \\v\\ is \\2^{53} - 1\\.

## Details

The Miller-Rabin primality test is a probabilistic algorithm that makes
heavy use of [modular
exponentiation](https://en.wikipedia.org/wiki/Modular_exponentiation).
At the heart of modular exponentiation is the ability to accurately
obtain the remainder of the product of two numbers \\\pmod p\\.

With the gmp library, producing accurate calculations for problems like
this is trivial because of the nature of the multiple precision data
type. However, standard C++ does not afford this luxury and simply
relying on a strict translation would have limited this algorithm to
numbers less than \\\sqrt 2^{63} - 1\\ (N.B. We are taking advantage of
the signed 64-bit fixed width integer from the stdint library in C++. If
we were confined to base R, the limit would have been \\\sqrt 2^{53} -
1\\). RcppAlgos::isPrimeRcpp gets around this limitation with a [divide
and conquer](https://en.wikipedia.org/wiki/Divide_and_conquer_algorithm)
approach taking advantage of properties of arithmetic.

The problem we are trying to solve can be summarized as follows:

\$\$(x_1 \* x_2) \pmod p\$\$

Now, we rewrite \\x_2\\ as \\x_2 = y_1 + y_2 + \dots + y_n\\, so that we
obtain:

\$\$(x_1 \* y_1) \pmod p + (x_1 \* y_2) \pmod p + \dots + (x_1 \* y_n)
\pmod p\$\$

Where each product \\(x_1 \* y_j)\\ for \\j \<= n\\ is smaller than the
original \\x_1 \* x_2\\. With this approach, we are now capable of
handling much larger numbers. Many details have been omitted for
clarity.

For a more in depth examination of this topic see [Accurate Modular
Arithmetic with Double
Precision](https://codereview.stackexchange.com/questions/186751/accurate-modular-arithmetic-with-double-precision).

## Value

Returns a named/unnamed logical vector. If an index is `TRUE`, the
number at that index is prime, otherwise the number is composite.

## References

- [THE MILLER-RABIN
  TEST](https://www.math.uconn.edu/~kconrad/blurbs/ugradnumthy/millerrabin.pdf)

  - Conrad, Keith. "THE MILLER-RABIN TEST."
    https://www.math.uconn.edu/~kconrad/blurbs/ugradnumthy/millerrabin.pdf.

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## See also

[`primeFactorize`](https://jwood000.github.io/RcppAlgos/reference/primeFactorize.md),
[`isprime`](https://rdrr.io/pkg/gmp/man/isprime.html)

## Examples

``` r
## check the primality of a single number
isPrimeRcpp(100)
#> [1] FALSE

## check the primality of every number in a vector
isPrimeRcpp(1:100)
#>   [1] FALSE  TRUE  TRUE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
#>  [13]  TRUE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
#>  [25] FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE
#>  [37]  TRUE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
#>  [49] FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
#>  [61]  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
#>  [73]  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
#>  [85] FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#>  [97]  TRUE FALSE FALSE FALSE

set.seed(42)
mySamp <- sample(10^13, 10)

## return named vector for easy identification
isPrimeRcpp(mySamp, namedVector = TRUE)
#> 5053637821668 4945473353752  708576091668 3861127950937 5611435811813 
#>         FALSE         FALSE         FALSE          TRUE         FALSE 
#> 8651983869445 2062476240920 5651694855126 9421521179996 9639618136290 
#>         FALSE         FALSE         FALSE         FALSE         FALSE 

## Using nThreads
system.time(isPrimeRcpp(mySamp, nThreads = 2))
#>    user  system elapsed 
#>       0       0       0 
```
