# Count SemiPrimes With an Optional Partial-Sieve Cutoff

Count positive integers less than or equal to \\n\\ with exactly two
prime factors, counted with multiplicity, excluding the first \\a\\
primes as factors. This is the correction term \\P_2(n, a)\\ used in
Meissel and Lehmer prime counting.

## Usage

``` r
semiPrimeCount(n, a = 0, nThreads = NULL)
```

## Arguments

- n:

  A single nonnegative whole number, at most \\2^{53} - 1\\.

- a:

  A single nonnegative whole number, at most \\2^{53} - 1\\. This is the
  number of initial primes excluded as factors, not a prime value. The
  default, `a = 0`, permits every prime. For example, `a = 1` excludes
  the prime 2, so only odd semiprimes are counted.

- nThreads:

  Requested number of threads. The default is `NULL` for serial
  execution. Requests are capped by the available hardware and the
  useful work in the calculation. Small or simple cases may use fewer
  threads. This argument is validated even when the result is immediate.

## Details

Logical inputs `FALSE` and `TRUE` are converted to 0 and 1,
respectively. The thread count must be positive, so `nThreads = FALSE`
is not supported.

### Mathematical definition

Let \\p_i\\ denote the \\i\\th prime. The function counts products \\p_i
p_j \le n\\ with \\a \< i \le j\\. The ordering counts each integer
once: \\2 \cdot 3\\ and \\3 \cdot 2\\ do not give two different
semiprimes. Equal factors are allowed, so prime squares such as 4, 9,
and 25 are included when their prime factor is permitted.

When \\a = 0\\, the excluded set is empty; no zeroth prime is needed.
The ten semiprimes at most 30 are `4, 6, 9, 10, 14, 15, 21, 22, 25, 26`.
Excluding the first prime leaves `9, 15, 21, 25`, so
`semiPrimeCount(30, 1)` is 4.

Write \\\pi(t)\\ for the number of primes at most \\t\\, and set \\b =
\pi(\lfloor \sqrt n \rfloor)\\. For a fixed smaller factor \\p_i\\, the
larger factor has an index from \\i\\ through \\\pi(\lfloor n / p_i
\rfloor)\\. Consequently, \$\$P_2(n, a) = \sum\_{i=a+1}^{b}
\left\[\pi(\lfloor n / p_i \rfloor) - i + 1\right\].\$\$ A sum with no
terms is zero. In particular, the answer is zero for \\n \< 4\\ or \\a
\ge b\\. Equivalently, the smallest permitted prime square must be at
most \\n\\ for any semiprime to qualify.

### Relation to prime counting

For Meissel's cutoff \\a = \pi(\lfloor n^{1/3} \rfloor)\\, the partial
sieve leaves 1, the primes beyond the cutoff, and qualifying semiprimes.
Thus, for \\n \ge 1\\, \$\$\pi(n) = \Phi(n, a) + a - 1 - P_2(n, a).\$\$
The function
[`primePhi`](https://jwood000.github.io/RcppAlgos/reference/primePhi.md)
supplies \\\Phi(n, a)\\. Lehmer's smaller cutoff also requires
[`triPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/triPrimeCount.md);
see
[`primeCount`](https://jwood000.github.io/RcppAlgos/reference/primeCount.md)
for the full formula.

## Note

The calculation prepares prime and prime-count tables with storage
proportional to \\\sqrt n\\. Larger prime-count queries in the sum use
the shared partial-sieve evaluator or a segmented sieve, selected
according to the input, cutoff, and useful thread budget. The segmented
sieve uses bounded buffers for each worker; the function does not
allocate a dense table through \\n / p\_{a+1}\\.

At the largest supported \\n\\, the square-root table has limit
94,906,265. Its dense prime-count entries alone require \\4(L + 1)\\
bytes for limit \\L\\, about 362 MiB at the maximum, excluding prime
vectors, the shared Phi cache, temporary storage, and allocation
overhead. Workers within one call share the prepared tables; separate
calls prepare their own tables. A valid numeric input may still require
substantial time or memory.

## Value

A single whole number giving \\P_2(n, a)\\. The result has type
`integer` when it is at most `.Machine$integer.max`, and type `double`
otherwise. The count is exact throughout the supported range.

## Author

Joseph Wood

## See also

[`primePhi`](https://jwood000.github.io/RcppAlgos/reference/primePhi.md),
[`triPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/triPrimeCount.md),
[`primeCount`](https://jwood000.github.io/RcppAlgos/reference/primeCount.md),
[`primeSieve`](https://jwood000.github.io/RcppAlgos/reference/primeSieve.md)

## Examples

``` r
## Count all semiprimes, then exclude 2 as a factor
semiPrimeCount(30)    ## 10
#> [1] 10
semiPrimeCount(30, 1) ## 4
#> [1] 4

## Equal prime factors are included
semiPrimeCount(4)     ## 1, from 2 * 2
#> [1] 1
semiPrimeCount(3)     ## 0
#> [1] 0
semiPrimeCount(0)     ## 0
#> [1] 0
semiPrimeCount(30, 100) ## 0
#> [1] 0

## A modest calculation with a requested thread budget
semiPrimeCount(10^6, nThreads = 2)
#> [1] 210035
```
