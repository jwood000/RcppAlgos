# Find the nth Prime Number

Find the prime number at a given position in the increasing sequence of
primes. Positions start at 1, so the first prime is 2.

## Usage

``` r
nthPrime(n, nThreads = NULL)
```

## Arguments

- n:

  A single positive whole number, at most 252,252,704,148,404. This is
  the position of the prime to find, not an upper bound on the prime's
  value.

- nThreads:

  Requested number of threads. The default is `NULL` for serial
  execution. Requests are capped by the available hardware and the
  useful work in the calculation. Small or simple cases may use fewer
  threads. This argument is validated even when the result is immediate.

## Details

Let \\p_n\\ denote the \\n\\th prime. The result is the unique integer
\\p_n\\ satisfying \$\$\pi(p_n - 1) = n - 1, \qquad \pi(p_n) = n,\$\$
where \\\pi(x)\\ counts the primes less than or equal to \\x\\. For
example, `nthPrime(25)` returns 97, whereas `primeCount(25)` returns 9.

Small positions use a lookup table. Otherwise, the calculation starts
with an approximation to \\p_n\\, obtains an exact prime count at that
location, and finishes with a segmented sieve. At most one additional
whole-range prime count is used. Counts over nearby intervals may then
refine the search location without recounting all preceding primes.
Approximation errors affect the amount of work, not the correctness of
the result: the sieve continues until it finds the required prime.

The whole-range counts use Legendre, Meissel, or Lehmer, with the same
size and thread thresholds used to choose among those methods in
[`primeCount`](https://jwood000.github.io/RcppAlgos/reference/primeCount.md).
Prepared prime and prime-count tables are reused for these counts, then
released before interval sieving. Simple size bands, based on the search
location rather than the requested position, determine when to switch to
sieving. Interval counting and final selection reuse the same sieve as
[`primeSieve`](https://jwood000.github.io/RcppAlgos/reference/primeSieve.md).
Requested threads apply to both counting and sieving, subject to the
useful work available.

The largest supported position is \\\pi(2^{53} - 1)\\, so the returned
prime is at most \\2^{53} - 1\\ and is represented exactly. The table
memory grows with the square root of the search location. Near the upper
limit, the dense prime-count table alone occupies about 362 MiB, in
addition to the prime vector, shared Phi cache, and worker storage.
Final selection retains generated primes from a window of at most one
million integer positions at a time. Interval counting retains no
generated prime values; both sieve stages also require base-prime and
worker storage.

Logical inputs `TRUE` are converted to 1. Missing values, `FALSE`,
non-finite values, fractional values, and vectors of length other than 1
are not supported. `nThreads = NULL` is supported as described above.

## Value

A single prime number. The result is of type `integer` when it is at
most `.Machine$integer.max`, and `numeric` otherwise.

## See also

[`primeCount`](https://jwood000.github.io/RcppAlgos/reference/primeCount.md),
[`primeSieve`](https://jwood000.github.io/RcppAlgos/reference/primeSieve.md),
[`isPrimeRcpp`](https://jwood000.github.io/RcppAlgos/reference/isPrimeRcpp.md)

## References

Axler, C. (2018). *New estimates for the nth prime number*.
<https://arxiv.org/abs/1706.03651>

Walisch, K. *primecount*: prime counting and nth-prime algorithms.
<https://github.com/kimwalisch/primecount>

## Examples

``` r
nthPrime(1)
#> [1] 2
nthPrime(25)
#> [1] 97
nthPrime(1000000)
#> [1] 15485863

## The result has exactly n primes at or below it
n <- 1000
p <- nthPrime(n)
primeCount(p)
#> [1] 1000
primeCount(p - 1)
#> [1] 999

## Use multiple threads for counting and sieving when useful
nthPrime(1000000, nThreads = 2)
#> [1] 15485863
```
