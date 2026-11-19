# Prime Counting Function $`\pi(x)`$

Count the prime numbers less than or equal to $`n`$, using automatic
selection or an explicitly selected Legendre, Meissel, or Lehmer method.

## Usage

``` r
primeCount(n, nThreads = NULL, method = "auto")
```

## Arguments

- n:

  A single positive whole number, at most $`2^{53} - 1`$.

- nThreads:

  Requested number of threads. The default is `NULL` for serial
  execution. Requests are capped by the available hardware and the
  useful work in the calculation. Small or simple cases may use fewer
  threads. For `n < 100000`, this argument is ignored and execution is
  serial.

- method:

  One of `"auto"`, `"legendre"`, `"meissel"`, or `"lehmer"`. Matching is
  case-sensitive and accepts unique abbreviations, as in
  [`match.arg`](https://rdrr.io/r/base/match.arg.html). `NULL` selects
  `"auto"`.

## Details

### Method selection

The current `"auto"` policy selects a method using `n` and the usable
thread budget, which is the requested number of threads capped by the
available hardware:

- For `n < 100000`, use small lookup values or a sieve.

- For `100000 <= n < 10^10`, use Legendre's formula.

- For `n >= 10^10`, use Meissel's formula, except that inputs
  `n >= 10^12` with at least six usable threads use Lehmer's formula.

These choices are based on first-call and warmed-call benchmarks. They
do not depend on the current cache contents, and they do not guarantee
the fastest method for every input or machine. The policy may change
with further benchmarking. In particular, a method that improves
repeated calls can have a slower first call.

An explicit method uses the named algorithm for every accepted `n`,
including small inputs. The `method` argument follows `nThreads`, so
existing positional calls such as `primeCount(n, 2)` retain their
meaning. All methods use serial execution when `nThreads = NULL`.

### Counting formulas

Let $`p_a`$ be the $`a`$th prime, and let $`\Phi(n, a)`$ count the
positive integers at most $`n`$ that are not divisible by any of the
first $`a`$ primes. Its recurrence is
``` math
\Phi(n, a) = \Phi(n, a - 1) - \Phi(\lfloor n / p_a \rfloor, a - 1)
```
with $`\Phi(n, 0) = n`$.

For $`k = 2`$ or $`3`$, let $`P_k(n, a)`$ count integers at most $`n`$
with exactly $`k`$ prime factors, none among the first $`a`$ primes.
Factors are counted with multiplicity: $`p^2`$ is included in $`P_2`$,
and $`p^3`$ and $`p^2 q`$ are included in $`P_3`$ when their factors
exceed the cutoff. Nondecreasing factor order counts each integer once.
When $`a = 0`$, all primes are allowed.

The methods use different cutoffs and correction terms:

- Legendre uses $`a = \pi(\lfloor \sqrt n \rfloor)`$:
  ``` math
  \pi(n) = \Phi(n, a) + a - 1.
  ```

- Meissel uses $`a = \pi(\lfloor n^{1/3} \rfloor)`$:
  ``` math
  \pi(n) = \Phi(n, a) + a - 1 - P_2(n, a).
  ```

- Lehmer uses $`a = \pi(\lfloor n^{1/4} \rfloor)`$:
  ``` math
  \pi(n) = \Phi(n, a) + a - 1 - P_2(n, a) - P_3(n, a).
  ```

The Legendre implementation and shared $`\Phi`$ calculation build on Kim
Walisch's work in
[primecount](https://github.com/kimwalisch/primecount), including cached
partial-sieve values, small-prime periodic tables, and prime-count
lookup shortcuts.

## Note

The counting methods prepare tables with storage proportional to
$`\sqrt n`$. At the maximum input $`2^{53} - 1`$, the dense prime-count
table requires about 362 MiB and the prime-vector reservation requests
about 42 MiB, before caches and other overhead. These are calculated
storage estimates, not measured peak memory. Workers within one call
share the tables; concurrent calls prepare separate tables. The
supported numeric range does not guarantee that every input is practical
for available memory.

## References

- Tomás Oliveira e Silva, [Computing $`\pi(x)`$: the combinatorial
  method](https://sweet.ua.pt/tos/bib/5.4.pdf), Revista do DETUA, vol.
  4, no. 6, March 2006, p. 761.

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## Value

A single whole number giving the number of primes less than or equal to
$`n`$. The result has type `integer` when it is at most
`.Machine$integer.max`, and type `double` otherwise. The count is exact
throughout the supported range.

## Author

Joseph Wood

## See also

[`primeSieve`](https://jwood000.github.io/RcppAlgos/reference/primeSieve.md),
[`nthPrime`](https://jwood000.github.io/RcppAlgos/reference/nthPrime.md),
[`primePhi`](https://jwood000.github.io/RcppAlgos/reference/primePhi.md),
[`semiPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/semiPrimeCount.md),
[`triPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/triPrimeCount.md)

## Examples

``` r
## Get the number of primes less than or equal to a billion
primeCount(10^9)
#> [1] 50847534

## Select a method explicitly
primeCount(10^6, method = "legendre")
#> [1] 78498
primeCount(10^6, method = "meissel")
#> [1] 78498
primeCount(10^6, method = "lehmer")
#> [1] 78498

## Using nThreads
system.time(primeCount(10^10, nThreads = 2, method = "lehmer"))
#>    user  system elapsed 
#>   0.011   0.000   0.006 
```
