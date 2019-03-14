\name{eulerPhiSieve}
\alias{eulerPhiSieve}
\title{
Apply Euler's Phi Function to Every Element in a Range
}
\description{
Sieve that generates the number of coprime elements for every number between \code{bound1} and \code{bound2} (if supplied) or all numbers up to \code{bound1}. This is equivalent to applying Euler's phi function (often written as \eqn{\phi(x)}) to every number in a given range.
}
\usage{
eulerPhiSieve(bound1, bound2 = NULL, namedVector = FALSE, nThreads = NULL)
}
\arguments{
    \item{bound1}{Positive integer or numeric value.}
    \item{bound2}{Positive integer or numeric value.}
    \item{namedVector}{Logical flag. If \code{TRUE}, a named vector is returned. The default is \code{FALSE}.}
    \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
}
\details{
For the simple case (i.e. when \code{bound2 = NULL}), this algorithm first generates all primes up to \eqn{n} via the sieve of Eratosthenes. We use these primes to sieve over the sequence \code{1:n}, dividing each value by \eqn{p}, creating a temporary value that will be subtracted from the original value at each index (i.e. equivalent to multiply each index by \eqn{(1 - 1/p)} but more efficient as we don't have to deal with floating point numbers). The case when \code{is.null(bound2) = FALSE} is more complicated but the basic ideas still hold.

This function is very useful when you need to calculate Euler's phi function for many numbers in a range as performing this calculation on the fly can be computationally expensive.

This algorithm benefits greatly from the fast integer division library 'libdivide'. The following is from \href{http://libdivide.com/}{http://libdivide.com/}:
\itemize{\item{\dQuote{\emph{libdivide allows you to replace expensive integer divides with comparatively cheap multiplication and bitshifts. Compilers usually do this, but only when the divisor is known at compile time. libdivide allows you to take advantage of it at runtime. The result is that integer division can become faster - a lot faster.}}}}
}
\value{
Returns a named/unnamed integer vector if \code{max(bound1, bound2)} \eqn{< 2^{31}}{< 2^31}, or a numeric vector otherwise.
}
\author{
Joseph Wood
}
\note{
The maximum allowed value is \eqn{2^{53} - 1}{2^53 - 1}.
}
\references{
\itemize{
\item{\href{https://en.wikipedia.org/wiki/Euler\%27s_totient_function}{Euler's totient function}}
\item{\href{http://ridiculousfish.com/}{ridiculousfish (author of libdivide)}}
\item{\href{https://github.com/ridiculousfish/libdivide}{github.com/ridiculousfish/libdivide}}
\item{\href{https://en.wikipedia.org/wiki/Double-precision_floating-point_format}{53-bit significand precision}}
}
}
\seealso{
\code{\link[numbers]{eulersPhi}}
}
\examples{
## Generate some random data
set.seed(496)
mySamp <- sample(10^6, 5*10^5)

## Generate number of coprime elements for many numbers
system.time(myPhis <- eulerPhiSieve(10^6))

## Now use result in algorithm
for (s in mySamp) {
    sPhi <- myPhis[s]
    ## Continue algorithm
}

## See https://projecteuler.net
system.time(which.max((1:10^6)/eulerPhiSieve(10^6)))

## Generating number of coprime elements
## for every number in a range is no problem
system.time(myPhiRange <- eulerPhiSieve(10^13, 10^13 + 10^6))

## Returning a named vector
eulerPhiSieve(10, 20, namedVector = TRUE)
eulerPhiSieve(10, namedVector = TRUE)

## Using nThreads
system.time(eulerPhiSieve(1e5, 2e5, nThreads = 2))
}