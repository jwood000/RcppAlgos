# Generate Prime Numbers

Implementation of the segmented sieve of Eratosthenes with wheel
factorization. Generates all prime numbers between `bound1` and `bound2`
(if supplied) or all primes up to `bound1`. See this stackoverflow post
for an analysis on prime number generation efficiency in R: [Generate a
list of primes up to a certain
number](https://stackoverflow.com/a/48313378/4408538)

The fundamental concepts of this algorithm are based off of the
implementation by Kim Walisch found here:
[kimwalisch/primesieve](https://github.com/kimwalisch/primesieve).

## Usage

``` r
primeSieve(bound1, bound2 = NULL, nThreads = NULL)
```

## Arguments

- bound1:

  Positive integer or numeric value.

- bound2:

  Positive integer or numeric value.

- nThreads:

  Specific number of threads to be used. The default is `NULL`.

## Note

- It does not matter which bound is larger as the resulting primes will
  be between `min(bound1, bound2)` and `max(bound1, bound2)` if `bound2`
  is provided.

- The maximum value for either of the bounds is \\2^{53} - 1\\.

## Details

At the heart of this algorithm is the traditional sieve of Eratosthenes
(i.e. given a [prime](https://en.wikipedia.org/wiki/Prime_number) \\p\\,
mark all multiples of \\p\\ as
[composite](https://en.wikipedia.org/wiki/Composite_number)), however
instead of sieving the entire interval, we only consider small
sub-intervals. The benefits of this method are two fold:

1.  Reduction of the [space
    complexity](https://en.wikipedia.org/wiki/DSPACE) from \\O(n)\\, for
    the traditional sieve, to \\O(\sqrt n)\\

2.  Reduction of [cache
    misses](https://en.wikipedia.org/wiki/CPU_cache#Cache_miss)

The latter is of particular importance as cache memory is much more
efficient and closer in proximity to the CPU than [main
memory](https://en.wikipedia.org/wiki/Computer_data_storage#Primary_storage).
Reducing the size of the sieving interval allows for more effective
utilization of the cache, which greatly impacts the overall efficiency.

Another optimization over the traditional sieve is the utilization of
wheel factorization. With the traditional sieve of Eratosthenes, you
typically check every odd index of your logical vector and if the value
is true, you have found a prime. With wheel factorization using the
first four primes (i.e. 2, 3, 5, and 7) to construct your wheel (i.e.
210 wheel), you only have to check indices of your logical vector that
are coprime to 210 (i.e. the product of the first four primes). As an
example, with \\n = 10000\\ and a 210 wheel, you only have to check 2285
indices vs. 5000 with the classical implementation.

## Value

Returns an integer vector if `max(bound1, bound2)` \\\< 2^{31}\\, or a
numeric vector otherwise.

## References

- [primesieve (Fast C/C++ prime number
  generator)](https://github.com/kimwalisch/primesieve)

- [Sieve of
  Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes)

- [Wheel
  factorization](https://en.wikipedia.org/wiki/Wheel_factorization)

- [53-bit significand
  precision](https://en.wikipedia.org/wiki/Double-precision_floating-point_format)

## Author

Joseph Wood

## Examples

``` r
## Primes up to a thousand
primeSieve(100)
#>  [1]  2  3  5  7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97

## Primes between 42 and 17
primeSieve(42, 17)
#> [1] 17 19 23 29 31 37 41

## Equivalent to
primeSieve(17, 42)
#> [1] 17 19 23 29 31 37 41

## Primes up to one hundred million in no time
system.time(primeSieve(10^8))
#>    user  system elapsed 
#>   0.111   0.005   0.115 

## options(scipen = 50)
## Generate large primes over interval
system.time(myPs <- primeSieve(10^13+10^6, 10^13))
#>    user  system elapsed 
#>    0.01    0.00    0.01 
## Object created is small
object.size(myPs)
#> 267696 bytes

## Using nThreads
system.time(primeSieve(1e7, nThreads = 2))
#>    user  system elapsed 
#>   0.012   0.000   0.007 
```
