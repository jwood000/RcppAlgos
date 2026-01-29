# Prime Counting Function \\\pi(x)\\

[Prime counting
function](https://en.wikipedia.org/wiki/Prime-counting_function) for
counting the prime numbers less than an integer, \\n\\, using Legendre's
formula. It is based on the the algorithm developed by Kim Walisch found
here: [kimwalisch/primecount](https://github.com/kimwalisch/primecount).

## Usage

``` r
primeCount(n, nThreads = NULL)
```

## Arguments

- n:

  Positive number

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

## Details

[Legendre's
Formula](https://mathworld.wolfram.com/LegendresFormula.html) for
counting the number of primes less than \\n\\ makes use of the
[inclusion-exclusion
principle](https://en.wikipedia.org/wiki/Inclusion-exclusion_principle)
to avoid explicitly counting every prime up to \\n\\. It is given by:
\$\$\pi(x) = \pi(\sqrt x) + \Phi(x, \sqrt x) - 1\$\$ Where \\\Phi(x,
a)\\ is the number of positive integers less than or equal to \\x\\ that
are relatively prime to the first \\a\\ primes (i.e. not divisible by
any of the first \\a\\ primes). It is given by the recurrence relation
(\\p_a\\ is the \\ath\\ prime (e.g. \\p_4 = 7\\)): \$\$\Phi(x, a) =
\Phi(x, a - 1) + \Phi(x / p_a, a - 1)\$\$ This algorithm implements five
modifications developed by Kim Walisch for calculating \\\Phi(x, a)\\
efficiently.

1.  Cache results of \\\Phi(x, a)\\

2.  Calculate \\\Phi(x, a)\\ using \\\Phi(x, a) = (x / pp) \* \phi(pp) +
    \Phi(x mod pp, a)\\ if \\a \<= 6\\

    - \\pp = 2 \* 3 \* ... \* \\ `prime[a]`

    - \\\phi(pp) = (2 - 1) \* (3 - 1) \* ... \* \\ \\(\\`prime[a]` \\-
      1)\\ (i.e. Euler's totient function)

3.  Calculate \\\Phi(x, a)\\ using \\\pi(x)\\ lookup table

4.  Calculate all \\\Phi(x, a) = 1\\ upfront

5.  Stop recursion at \\6\\ if \\\sqrt x \>= 13\\ or \\\pi(\sqrt x)\\
    instead of \\1\\

## Note

The maximum value of \\n\\ is \\2^{53} - 1\\

## References

- [Computing \\\pi(x)\\: the combinatorial
  method](https://sweet.ua.pt/tos/bib/5.4.pdf)

  - Tomás Oliveira e Silva, Computing pi(x): the combinatorial method,
    Revista do DETUA, vol. 4, no. 6, March 2006, p. 761.
    https://sweet.ua.pt/tos/bib/5.4.pdf

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## Value

Whole number representing the number of prime numbers less than or equal
to \\n\\.

## Author

Joseph Wood

## See also

[`primeSieve`](https://jwood000.github.io/RcppAlgos/reference/primeSieve.md)

## Examples

``` r
## Get the number of primes less than a billion
primeCount(10^9)
#> [1] 50847534

## Using nThreads
system.time(primeCount(10^10, nThreads = 2))
#>    user  system elapsed 
#>   0.019   0.002   0.015 
```
