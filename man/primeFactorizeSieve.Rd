\name{primeFactorizeSieve}
\alias{primeFactorizeSieve}
\title{
Generate Prime Factorization for Numbers in a Range
}
\description{
Generates the prime factorization of all numbers between \code{bound1} and \code{bound2} (if supplied) or all numbers up to \code{bound1}.
}
\usage{
primeFactorizeSieve(bound1, bound2 = NULL, namedList = FALSE, nThreads = NULL)
}
\arguments{
    \item{bound1}{Positive integer or numeric value.}
    \item{bound2}{Positive integer or numeric value.}
    \item{namedList}{Logical flag. If \code{TRUE}, a named list is returned. The default is \code{FALSE}.}
    \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
}
\details{
This function is useful when many prime factorizations are needed. Instead of generating the prime factorization on the fly, one can reference the indices/names of the generated list.

This algorithm benefits greatly from the fast integer division library 'libdivide'. The following is from \href{http://libdivide.com/}{http://libdivide.com/}:
\itemize{\item{\dQuote{\emph{libdivide allows you to replace expensive integer divides with comparatively cheap multiplication and bitshifts. Compilers usually do this, but only when the divisor is known at compile time. libdivide allows you to take advantage of it at runtime. The result is that integer division can become faster - a lot faster.}}}}
}
\value{
Returns a named/unnamed list of integer vectors if \code{max(bound1, bound2)} \eqn{< 2^{31}}{< 2^31}, or a list of numeric vectors otherwise.
}
\note{
The maximum value for either of the bounds is \eqn{2^{53} - 1}{2^53 - 1}.
}
\references{
\itemize{
\item{\href{https://en.wikipedia.org/wiki/Prime_factor}{Prime Factor}}
\item{\href{http://ridiculousfish.com/}{ridiculousfish (author of libdivide)}}
\item{\href{https://github.com/ridiculousfish/libdivide}{github.com/ridiculousfish/libdivide}}
\item{\href{https://en.wikipedia.org/wiki/Double-precision_floating-point_format}{53-bit significand precision}}
}
}
\author{
Joseph Wood
}
\seealso{
\code{\link{primeFactorize}}, \code{\link{divisorsSieve}}, \code{\link[gmp]{factorize}}, \code{\link[numbers]{primeFactors}}
}
\examples{
## Generate some random data
set.seed(28)
mySamp <- sample(10^5, 5*10^4)

## Generate prime factorizations up 
## to 10^5 (max element from mySamp)
system.time(allPFacs <- primeFactorizeSieve(10^5))

## Use generated prime factorization for further
## analysis by accessing the index of allPFacs
for (s in mySamp) {
    pFac <- allPFacs[[s]]
    ## Continue algorithm
}

## Generating prime factorizations over
## a range is efficient as well
system.time(primeFactorizeSieve(10^12, 10^12 + 10^5))

## Set 'namedList' to TRUE to return a named list
primeFactorizeSieve(27, 30, namedList = TRUE)

## Using nThreads
system.time(primeFactorizeSieve(1e4, 5e4, nThreads = 2))
}
