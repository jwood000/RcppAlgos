\name{numDivisorSieve}
\alias{numDivisorSieve}
\title{
Apply Divisor Function to Every Element in a Range
}
\description{
Sieve that generates the number of divisors for every number between \code{bound1} and \code{bound2} (if supplied) or all numbers up to \code{bound1}. This is equivalent to applying the divisor function (often written as \eqn{\sigma(x)}) to every number in a given range.
}
\usage{
numDivisorSieve(bound1, bound2 = NULL, namedVector = FALSE, nThreads = NULL)
}
\arguments{
    \item{bound1}{Positive integer or numeric value.}
    \item{bound2}{Positive integer or numeric value.}
    \item{namedVector}{Logical flag. If \code{TRUE}, a named vector is returned. The default is \code{FALSE}.}
    \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
}
\details{
Simple and efficient sieve that calculates the number of divisors for every number in a given range. This function is very useful when you need to calculate the number of divisors for many numbers.

This algorithm benefits greatly from the fast integer division library 'libdivide'. The following is from \href{http://libdivide.com/}{http://libdivide.com/}:
\itemize{\item{\dQuote{\emph{libdivide allows you to replace expensive integer divides with comparatively cheap multiplication and bitshifts. Compilers usually do this, but only when the divisor is known at compile time. libdivide allows you to take advantage of it at runtime. The result is that integer division can become faster - a lot faster.}}}}
}
\value{
Returns a named/unnamed integer vector
}
\author{
Joseph Wood
}
\note{
The maximum allowed value is \eqn{2^{53} - 1}{2^53 - 1}.
}
\references{
\itemize{
\item{\href{https://en.wikipedia.org/wiki/Divisor_function}{Divisor function}}
\item{\href{http://ridiculousfish.com/}{ridiculousfish (author of libdivide)}}
\item{\href{https://github.com/ridiculousfish/libdivide}{github.com/ridiculousfish/libdivide}}
\item{\href{https://en.wikipedia.org/wiki/Double-precision_floating-point_format}{53-bit significand precision}}}
}
\examples{
## Generate some random data
set.seed(8128)
mySamp <- sample(10^6, 5*10^5)

## Generate number of divisors for
## every number less than a million
system.time(mySigmas <- numDivisorSieve(10^6))

## Now use result in algorithm
for (s in mySamp) {
    sSig <- mySigmas[s]
    ## Continue algorithm
}

## Generating number of divisors for every
## number in a range is no problem
system.time(sigmaRange <- numDivisorSieve(10^13, 10^13 + 10^6))

## Returning a named vector
numDivisorSieve(10, 20, namedVector = TRUE)
numDivisorSieve(10, namedVector = TRUE)

## Using nThreads
system.time(numDivisorSieve(1e5, 2e5, nThreads = 2))
}