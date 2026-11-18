# Count Integers Remaining After a Partial Prime Sieve

Count the positive integers less than or equal to \\n\\ that are not
divisible by any of the first \\a\\ primes. This is the partial-sieve
function \\\Phi(n, a)\\ used in combinatorial prime counting, not
Euler's totient function.

## Usage

``` r
primePhi(n, a, nThreads = NULL)
```

## Arguments

- n:

  A single nonnegative whole number, at most \\2^{53} - 1\\.

- a:

  A single nonnegative whole number, at most \\2^{53} - 1\\. This is the
  number of initial primes to exclude, not a prime value or an upper
  bound on the excluded primes. For example, `a = 3` excludes multiples
  of 2, 3, and 5. With `a = 0`, no primes are excluded.

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

Write \\p_1 = 2, p_2 = 3, p_3 = 5, \ldots\\ for the increasing sequence
of primes. Then \\\Phi(n, a)\\ counts integers \\m\\ with \\1 \le m \le
n\\ for which none of \\p_1, \ldots, p_a\\ divides \\m\\. When \\a =
0\\, the excluded set is empty. In particular, \$\$\Phi(n, 0) = n,
\qquad \Phi(0, a) = 0.\$\$ For \\a \ge 1\\, removing multiples of the
next prime gives \$\$\Phi(n, a) = \Phi(n, a - 1) - \Phi(\lfloor n / p_a
\rfloor, a - 1).\$\$ The subtracted term counts the multiples of \\p_a\\
that survived removal of the preceding primes.

The integer 1 survives every partial sieve. Thus, when \\n \> 0\\, the
answer is at least 1, and it is exactly 1 when \\a \ge \pi(n)\\, where
\\\pi(n)\\ is the number of primes at most \\n\\. Composite numbers can
also survive: for example, `primePhi(30, 2)` includes 25 because neither
2 nor 3 divides it. This gives the ten survivors
`1, 5, 7, 11, 13, 17, 19, 23, 25, 29`.

Euler's totient \\\phi(n)\\ instead counts integers at most \\n\\ that
are coprime to \\n\\ itself. For that function, see
[`eulerPhiSieve`](https://jwood000.github.io/RcppAlgos/reference/eulerPhiSieve.md).

### Relation to prime counting

If \\a = \pi(\lfloor \sqrt n \rfloor)\\ and \\n \ge 1\\, every surviving
integer other than 1 is a prime larger than the excluded primes.
Legendre's formula is therefore \$\$\pi(n) = \Phi(n, a) + a - 1.\$\$
Smaller choices of \\a\\ leave composite survivors. Meissel and Lehmer
prime counting account for these using the correction terms exposed by
[`semiPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/semiPrimeCount.md)
and
[`triPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/triPrimeCount.md).
See
[`primeCount`](https://jwood000.github.io/RcppAlgos/reference/primeCount.md)
for the corresponding formulas.

## Note

The general calculation prepares prime and prime-count tables with
storage proportional to \\\sqrt n\\. Cases with \\a \le 6\\ use small
periodic tables instead of a full square-root table. Large values of
\\a\\ do not require preparing all of the first \\a\\ primes merely to
establish that the answer is 1.

The largest supported table limit is 94,906,265. Its dense prime-count
table alone requires \\4(L + 1)\\ bytes for limit \\L\\, about 362 MiB
at that maximum. Prime vectors, the shared Phi cache, temporary storage,
and allocation overhead require additional memory. Workers within one
call share the prepared tables; separate calls prepare their own tables.
The numeric input range does not guarantee that every calculation is
practical for the available resources.

## Value

A single whole number giving \\\Phi(n, a)\\. The result has type
`integer` when it is at most `.Machine$integer.max`, and type `double`
otherwise. The count is exact throughout the supported range.

## Author

Joseph Wood

## See also

[`primeCount`](https://jwood000.github.io/RcppAlgos/reference/primeCount.md),
[`semiPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/semiPrimeCount.md),
[`triPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/triPrimeCount.md),
[`eulerPhiSieve`](https://jwood000.github.io/RcppAlgos/reference/eulerPhiSieve.md),
[`primeSieve`](https://jwood000.github.io/RcppAlgos/reference/primeSieve.md)

## Examples

``` r
## Remove multiples of the first two primes, 2 and 3
primePhi(30, 2)  ## 10
#> [1] 10

## No excluded primes; an empty interval; all primes through n excluded
primePhi(30, 0)   ## 30
#> [1] 30
primePhi(0, 2)    ## 0
#> [1] 0
primePhi(30, 100) ## 1
#> [1] 1

## Recover the prime count using Legendre's formula
n <- 1000
a <- primeCount(floor(sqrt(n)))
primePhi(n, a) + a - 1 == primeCount(n)
#> [1] TRUE
```
