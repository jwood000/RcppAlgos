\name{eulerPhiSieve}
\alias{eulerPhiSieve}
\title{
Number of Coprime Integers for all Integers up to \eqn{n}
}
\description{
Rcpp sieve implementation that quickly generates the number of coprime integers up to \eqn{n}. This is equivalent to performing Euler's phi function on every number up to \eqn{n}.
}
\usage{
eulerPhiSieve(n = 100L)
}
\arguments{
  \item{n}{Positive integer or numeric value.}
}
\details{
This function first generates all primes up to \eqn{n} via the sieve of Eratosthenes. Subsequently we use these primes to sieve over the sequence \code{1:n}, multiplying each index by \eqn{(1 - 1/p)} for each \eqn{p}. Finally, we coerce the resulting elements to integers as the intermediate results are not integer types.

This function is very useful when you need to calculate Euler's phi function on many numbers as performing this calculation on the fly can be computationally expensive.
}
\value{
Returns an integer vector
}
\author{
Joseph Wood
}
\references{
\href{https://en.wikipedia.org/wiki/Euler\%27s_totient_function}{Euler's totient function}
}
\seealso{
\code{\link[numbers]{eulersPhi}}
}
\examples{
## Generate some random data
set.seed(496)
mySamp <- sample(10^6, 5*10^5)

## Quickly generate number of coprime elements for many numbers
system.time(myPhis <- eulerPhiSieve(10^6))

## Now use result in algorithm
for (s in mySamp) {
    sPhi <- myPhis[s]
    ## Continue algorithm
}

## See https://projecteuler.net
system.time(which.max((1:10^6)/eulerPhiSieve(10^6)))
}