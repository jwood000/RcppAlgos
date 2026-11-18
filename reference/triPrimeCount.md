# Count TriPrimes With an Optional Partial-Sieve Cutoff

Count positive integers less than or equal to \\n\\ with exactly three
prime factors, counted with multiplicity, excluding the first \\a\\
primes as factors. This is the correction term \\P_3(n, a)\\ used in
Lehmer prime counting. Its table requirements impose an additional
input-dependent limit described below.

## Usage

``` r
triPrimeCount(n, a = 0, nThreads = NULL)
```

## Arguments

- n:

  A single nonnegative whole number, at most \\2^{53} - 1\\. The
  additional table limit in ‘Details’ may restrict the accepted value
  for a given `a`. With the default `a = 0`, the largest accepted `n` is
  379,625,063.

- a:

  A single nonnegative whole number, at most \\2^{53} - 1\\. This is the
  number of initial primes excluded as factors, not a prime value. The
  default, `a = 0`, permits every prime. For example, `a = 1` excludes
  the prime 2, so only odd triprimes are counted. Larger cutoffs can
  permit larger `n` within the table limit.

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
p_j p_k \le n\\ with \\a \< i \le j \le k\\. Nondecreasing factor order
counts each integer once. Factors need not be distinct: prime cubes
\\p^3\\ and products \\p^2 q\\ are included, as are products of three
different primes.

When \\a = 0\\, the excluded set is empty; no zeroth prime is needed.
The seven triprimes at most 30 are `8, 12, 18, 20, 27, 28, 30`. For
example, 8 qualifies as \\2^3\\, while 16 does not qualify because
\\2^4\\ has four prime factors counted with multiplicity. Excluding 2
leaves only 27 among these seven values.

Write \\\pi(t)\\ for the number of primes at most \\t\\, and define
\$\$c = \pi(\lfloor n^{1/3} \rfloor), \qquad b_i = \pi(\lfloor \sqrt{n /
p_i} \rfloor).\$\$ The smallest factor has index at most \\c\\; for that
factor, the middle factor has index at most \\b_i\\. Once these two
factors are fixed, the largest factor has index from \\j\\ through
\\\pi(\lfloor n / (p_i p_j) \rfloor)\\. Hence, \$\$P_3(n, a) =
\sum\_{i=a+1}^{c} \sum\_{j=i}^{b_i} \left\[\pi(\lfloor n / (p_i p_j)
\rfloor) - j + 1\right\].\$\$ A sum with no terms is zero. In
particular, the answer is zero for \\n \< 8\\ or \\a \ge c\\: the
smallest permitted prime cube then exceeds \\n\\.

### Relation to prime counting

At Lehmer's cutoff \\a = \pi(\lfloor n^{1/4} \rfloor)\\, the surviving
composites have two or three permitted prime factors. For \\n \ge 1\\,
\$\$\pi(n) = \Phi(n, a) + a - 1 - P_2(n, a) - P_3(n, a).\$\$ The other
terms are available through
[`primePhi`](https://jwood000.github.io/RcppAlgos/reference/primePhi.md)
and
[`semiPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/semiPrimeCount.md).
See
[`primeCount`](https://jwood000.github.io/RcppAlgos/reference/primeCount.md)
for the complete counting methods.

### Table limit and memory use

The implementation first prepares tables through \\\lfloor n^{1/3}
\rfloor\\ to determine whether any qualifying triples exist. If they do,
let \\p\_{a+1}\\ be the smallest permitted prime. Direct prime-count
lookups require a table through \$\$q\_{\max} = \lfloor n / p\_{a+1}^2
\rfloor.\$\$ Every queried quotient is at most this value. The largest
supported table limit is 94,906,265; a nonempty calculation requiring a
larger table raises an error before that larger table is allocated. An
empty count does not require this second table.

With `a = 0`, the smallest permitted prime is 2, so the condition is
\\\lfloor n / 4 \rfloor \le 94906265\\. This explains the default upper
limit \\379625063\\. Excluding more primes can reduce the required table
size and allow larger inputs, up to the common numeric limit \\2^{53} -
1\\ when the table condition is met.

Table storage is proportional to \\\max(n^{1/3}, n / p\_{a+1}^2)\\ for
nonempty calculations. The dense prime-count entries alone use \\4(L +
1)\\ bytes for limit \\L\\, about 362 MiB at the maximum. This excludes
prime vectors, temporary tables, worker storage, allocation overhead,
and any shared Phi cache retained from other calls; it is not a
peak-memory bound. Workers within one call share the prepared lookup
tables. Inputs below the enforced table limit can still be expensive in
time or memory.

## Value

A single whole number giving \\P_3(n, a)\\. The result has type
`integer` when it is at most `.Machine$integer.max`, and type `double`
otherwise. The count is exact for accepted inputs.

## Author

Joseph Wood

## See also

[`primePhi`](https://jwood000.github.io/RcppAlgos/reference/primePhi.md),
[`semiPrimeCount`](https://jwood000.github.io/RcppAlgos/reference/semiPrimeCount.md),
[`primeCount`](https://jwood000.github.io/RcppAlgos/reference/primeCount.md),
[`primeSieve`](https://jwood000.github.io/RcppAlgos/reference/primeSieve.md)

## Examples

``` r
## Count all triprimes, then exclude 2 as a factor
triPrimeCount(30)    ## 7
#> [1] 7
triPrimeCount(30, 1) ## 1, from 3 * 3 * 3
#> [1] 1

## Repeated prime factors are included
triPrimeCount(8)      ## 1, from 2 * 2 * 2
#> [1] 1
triPrimeCount(7)      ## 0
#> [1] 0
triPrimeCount(0)      ## 0
#> [1] 0
triPrimeCount(30, 100) ## 0
#> [1] 0

## Keep exploration modest: the default table extends through floor(n / 4)
triPrimeCount(10^5, nThreads = 2)
#> [1] 25556
```
