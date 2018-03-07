\name{divisorsRcpp}
\alias{divisorsRcpp}
\title{
Vectorized Factorization (Complete)
}
\description{
Rcpp implementation that quickly generates the complete factorization for many numbers.
}
\usage{
divisorsRcpp(v = 100L, namedList = FALSE)
}
\arguments{
\item{v}{Vector of integers or numeric values. Non-integral values will be coerced to whole numbers.}
\item{namedList}{Logical flag. If \code{TRUE} and the \code{length(v) > 1}, a named list is returned. The default is \code{FALSE}.}
}
\details{
Highly efficient algorithm that builds on \code{\link{primeFactorize}} to quickly generate the complete factorization of many numbers.
}
\value{
\itemize{
\item{Returns an unnamed vector if \code{length(v) == 1} regardless of the value of \code{namedList}. If \eqn{v < 2^{31}}{v < 2^31}, the class of the returned vector will be integer, otherwise the class will be numeric.}
\item{If \code{length(v) > 1}, a named/unnamed list of vectors will be returned. If \code{max(bound1, bound2)} \eqn{< 2^{31}}{< 2^31}, the class of each vector will be integer, otherwise the class will be numeric.}
}
}
\references{
\itemize{
\item{\href{https://en.wikipedia.org/wiki/Divisor}{Divisor}}
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
\code{\link[numbers]{divisors}}, \code{\link{primeFactorize}}
}
\examples{
## Get the complete factorization of a single number
divisorsRcpp(10^8)

## Or get the complete factorization of many numbers
set.seed(29)
myVec <- sample(-1000000:1000000, 1000)
system.time(myFacs <- divisorsRcpp(myVec))

## Return named list
myFacsWithNames <- divisorsRcpp(myVec, namedList = TRUE)
}
