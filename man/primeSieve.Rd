\name{primeSieve}
\alias{primeSieve}
\title{
Generate or Count Prime Numbers
}
\description{
Implementation of the segmented sieve of Eratosthenes with wheel factorization. Generates all prime numbers between \code{bound1} and \code{bound2} (if supplied) or all primes up to \code{bound1}. Set \code{countOnly = TRUE} to return the number of primes in the same interval without storing the generated primes. Both endpoints are included. See this stackoverflow post for an analysis on prime number generation efficiency in R: \href{https://stackoverflow.com/a/48313378/4408538}{Generate a list of primes up to a certain number
}

The fundamental concepts of this algorithm are based off of the implementation by Kim Walisch found here: \href{https://github.com/kimwalisch/primesieve}{kimwalisch/primesieve}.
}
\usage{
primeSieve(bound1, bound2 = NULL, nThreads = NULL, countOnly = FALSE)
}
\arguments{
    \item{bound1}{Positive integer or numeric value.}
    \item{bound2}{Positive integer or numeric value.}
    \item{nThreads}{Requested number of threads. The default is
    \code{NULL} for serial execution. Requests are capped by the available
    hardware and the useful work in the calculation. Small or simple cases
    may use fewer threads.}
    \item{countOnly}{Logical value. If \code{TRUE}, returns only the number of
    primes in the interval. The default is \code{FALSE}, which returns the
    primes themselves.}
}
\note{
\itemize{
\item It does not matter which bound is larger as the resulting primes will be between \code{min(bound1, bound2)} and \code{max(bound1, bound2)} if \code{bound2} is provided.
\item For non-whole bounds, the lower bound is rounded up and the upper bound
is rounded down. These rules apply to both generation and counting.
\item The maximum value for either of the bounds is \eqn{2^{53} - 1}{2^53 - 1}.
}
}
\details{
When \code{countOnly = TRUE}, the same segmented sieve counts the surviving
primes directly. This avoids allocating a vector of all primes in the interval,
but the sieving primes, sieve segments, and other working data still require
memory. Parallel calculations may use additional working memory. Counting is
particularly useful for short intervals with large bounds.

Counting still uses the sieve when \code{bound2} is omitted. For the number of
primes up to a large value, use \code{\link{primeCount}}, which also provides
combinatorial prime-counting methods. Neither function counts 1 as prime.

At the heart of this algorithm is the traditional sieve of Eratosthenes (i.e. given a \href{https://en.wikipedia.org/wiki/Prime_number}{prime} \eqn{p}, mark all multiples of \eqn{p} as \href{https://en.wikipedia.org/wiki/Composite_number}{composite}), however instead of sieving the entire interval, we only consider small sub-intervals. The benefits of this method are two fold:
\enumerate{
\item Reduction of the \href{https://en.wikipedia.org/wiki/DSPACE}{space complexity} from \eqn{O(n)}, for the traditional sieve, to \eqn{O(\sqrt n)}
\item Reduction of \href{https://en.wikipedia.org/wiki/CPU_cache#Cache_miss}{cache misses}
}
The latter is of particular importance as cache memory is much more efficient and closer in proximity to the CPU than \href{https://en.wikipedia.org/wiki/Computer_data_storage#Primary_storage}{main memory}. Reducing the size of the sieving interval allows for more effective utilization of the cache, which greatly impacts the overall efficiency.

Another optimization over the traditional sieve is the utilization of wheel factorization. With the traditional sieve of Eratosthenes, you typically check every odd index of your logical vector and if the value is true, you have found a prime. With wheel factorization using the first four primes (i.e. 2, 3, 5, and 7) to construct your wheel (i.e. 210 wheel), you only have to check indices of your logical vector that are coprime to 210 (i.e. the product of the first four primes). As an example, with \eqn{n = 10000} and a 210 wheel, you only have to check 2285 indices vs. 5000 with the classical implementation.
}
\value{
When \code{countOnly = FALSE}, returns an integer vector if
\code{max(bound1, bound2)} \eqn{< 2^{31}}{< 2^31}, or a numeric vector otherwise.

When \code{countOnly = TRUE}, returns a single integer if the count is at most
\code{.Machine$integer.max}, or a single numeric value otherwise. The count is
exact throughout the supported range, and an interval containing no primes
returns \code{0L}. The count's type depends on the number of primes, not the
magnitude of the bounds.
}
\seealso{
\code{\link{primeCount}}, \code{\link{isPrimeRcpp}}
}
\references{
\itemize{
\item \href{https://github.com/kimwalisch/primesieve}{primesieve (Fast C/C++ prime number generator)}
\item \href{https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes}{Sieve of Eratosthenes}
\item \href{https://en.wikipedia.org/wiki/Wheel_factorization}{Wheel factorization}
\item \href{https://en.wikipedia.org/wiki/Double-precision_floating-point_format}{53-bit significand precision}
}
}
\author{
Joseph Wood
}
\examples{
## Primes up to a thousand
primeSieve(100)

## Primes between 42 and 17
primeSieve(42, 17)

## Equivalent to
primeSieve(17, 42)

## Count the primes without generating a vector
primeSieve(17, 42, countOnly = TRUE)

## Count primes in an interval with large bounds
primeSieve(1e10, 1e10 + 100, countOnly = TRUE)

## Count primes up to a hundred
primeSieve(100, countOnly = TRUE)

## Primes up to one hundred million in no time
system.time(primeSieve(10^8))

## options(scipen = 50)
## Generate large primes over interval
system.time(myPs <- primeSieve(10^13+10^6, 10^13))
## Object created is small
object.size(myPs)

## Using nThreads
system.time(primeSieve(1e7, nThreads = 2))
}
\keyword{prime}
\keyword{sieve}
\keyword{eratosthenes}
\keyword{number theory}
