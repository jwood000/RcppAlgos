\name{primeFactorize}
\alias{primeFactorize}
\title{
Vectorized Prime Factorization
}
\description{
Implementation of Pollard's rho algorithm for generating the prime factorization. The algorithm is based on the "factorize.c" source file from the gmp library found here \url{http://gmplib.org}.
}
\usage{
primeFactorize(v, namedList = FALSE, nThreads = NULL)
}
\arguments{
    \item{v}{Vector of integers or numeric values. Non-integral values will be cured to whole numbers.}
    \item{namedList}{Logical flag. If \code{TRUE} and the \code{length(v) > 1}, a named list is returned. The default is \code{FALSE}.}
    \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
}
\details{
As noted in the Description section above, this algorithm is based on the "factorize.c" source code from the gmp library. Much of the code in RcppAlgos::primeFactorize is a straightforward translation from multiple precision C data types to standard C++ data types. A crucial part of the algorithm's efficiency is based on quickly determining \href{https://en.wikipedia.org/wiki/Primality_test}{primality}, which is easily computed with gmp. However, with standard C++, this is quite challenging. Much of the research for RcppAlgos::primeFactorize was focused on developing an algorithm that could accurately and efficiently compute primality.

For more details, see the documentation for \code{\link{isPrimeRcpp}}.
}
\value{
\itemize{
\item{Returns an unnamed vector if \code{length(v) == 1} regardless of the value of \code{namedList}. If \eqn{v < 2^{31}}{v < 2^31}, the class of the returned vector will be integer, otherwise the class will be numeric.}
\item{If \code{length(v) > 1}, a named/unnamed list of vectors will be returned. If \code{max(bound1, bound2)} \eqn{< 2^{31}}{< 2^31}, the class of each vector will be integer, otherwise the class will be numeric.}
}
}
\references{
\itemize{
\item{\href{https://en.wikipedia.org/wiki/Pollard\%27s_rho_algorithm}{Pollard's rho algorithm}}
\item{\href{https://en.wikipedia.org/wiki/Miller-Rabin_primality_test}{Miller-Rabin primality test}}
\item{\href{https://codereview.stackexchange.com/questions/186751/accurate-modular-arithmetic-with-double-precision}{Accurate Modular Arithmetic with Double Precision}}
\item{\href{https://en.wikipedia.org/wiki/Double-precision_floating-point_format}{53-bit significand precision}}
}
}
\author{
Joseph Wood
}
\note{
The maximum value for each element in \eqn{v} is \eqn{2^{53} - 1}{2^53 - 1}.
}
\seealso{
\code{\link{primeFactorizeSieve}}, \code{\link[gmp]{factorize}}, \code{\link[numbers]{primeFactors}}
}
\examples{
## Get the prime factorization of a single number
primeFactorize(10^8)

## Or get the prime factorization of many numbers
set.seed(29)
myVec <- sample(-1000000:1000000, 1000)
system.time(pFacs <- primeFactorize(myVec))

## Return named list
pFacsWithNames <- primeFactorize(myVec, namedList = TRUE)

## Using nThreads
system.time(primeFactorize(myVec, nThreads = 2))
}
