\name{primeSieve}
\alias{primeSieve}
\title{
Generate Prime Numbers
}
\description{
Implementation of the segmented sieve of Eratosthenes with wheel factorization. Generates all prime numbers between \code{bound1} and \code{bound2} (if supplied) or all primes up to \code{bound1}. See this stackoverflow post for an analysis on prime number generation efficiency in R: \href{https://stackoverflow.com/a/48313378/4408538}{Generate a list of primes up to a certain number
}

The fundamental concepts of this algorithm are based off of the implementation by Kim Walisch found here: \href{https://github.com/kimwalisch/primesieve}{kimwalisch/primesieve}.
}
\usage{
primeSieve(bound1, bound2 = NULL, nThreads = NULL)
}
\arguments{
    \item{bound1}{Positive integer or numeric value.}
    \item{bound2}{Positive integer or numeric value.}
    \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
}
\note{
\itemize{
\item{It does not matter which bound is larger as the resulting primes will be between \code{min(bound1, bound2)} and \code{max(bound1, bound2)} if \code{bound2} is provided.}
\item{The maximum value for either of the bounds is \eqn{2^{53} - 1}{2^53 - 1}.}
}
}
\details{
At the heart of this algorithm is the traditional sieve of Eratosthenes (i.e. given a \href{https://en.wikipedia.org/wiki/Prime_number}{prime} \eqn{p}, mark all multiples of \eqn{p} as \href{https://en.wikipedia.org/wiki/Composite_number}{composite}), however instead of sieving the entire interval, we only consider small sub-intervals. The benefits of this method are two fold:
\enumerate{
\item Reduction of the \href{https://en.wikipedia.org/wiki/DSPACE}{space complexity} from \eqn{O(n)}, for the traditional sieve, to \eqn{O(\sqrt n)}
\item Reduction of \href{https://en.wikipedia.org/wiki/CPU_cache#Cache_miss}{cache misses}
}
The latter is of particular importance as cache memory is much more efficient and closer in proximity to the CPU than \href{https://en.wikipedia.org/wiki/Computer_data_storage#Primary_storage}{main memory}. Reducing the size of the sieving interval allows for more effective utilization of the cache, which greatly impacts the overall efficiency.

Another optimization over the traditional sieve is the utilization of wheel factorization. With the traditional sieve of Eratosthenes, you typically check every odd index of your logical vector and if the value is true, you have found a prime. With wheel factorization using the first four primes (i.e. 2, 3, 5, and 7) to construct your wheel (i.e. 210 wheel), you only have to check indices of your logical vector that are coprime to 210 (i.e. the product of the first four primes). As an example, with \eqn{n = 10000} and a 210 wheel, you only have to check 2285 indices vs. 5000 with the classical implementation.
}
\value{
Returns an integer vector if \code{max(bound1, bound2)} \eqn{< 2^{31}}{< 2^31}, or a numeric vector otherwise.
}
\references{
\itemize{
\item{\href{https://primesieve.org}{primesieve (Fast C/C++ prime number generator)}}
\item{\href{https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes}{Sieve of Eratosthenes}}
\item{\href{https://en.wikipedia.org/wiki/Wheel_factorization}{Wheel factorization}}
\item{\href{https://en.wikipedia.org/wiki/Double-precision_floating-point_format}{53-bit significand precision}}
}
}
\author{
Joseph Wood
}
\seealso{
\code{\link[numbers]{Primes}}
}
\examples{
## Primes up to a thousand
primeSieve(100)

## Primes between 42 and 17
primeSieve(42, 17)

## Equivalent to 
primeSieve(17, 42)

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