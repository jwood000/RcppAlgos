\name{numDivisorSieve}
\alias{numDivisorSieve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Number of Divisors for all Integers up to \eqn{n}
}
\description{
Rcpp sieve implementation that quickly generates the number of divisors for every integer up to \eqn{n}. This is equivalent to performing the divisor function (often written as \eqn{\sigma(n)}) for every number up to \eqn{n}.
}
\usage{
numDivisorSieve(n = 100L)
}
\arguments{
  \item{n}{Positive integer or numeric value.}
}
\details{
Very simple sieve that quickly calculates the number of divisors for every integer up to \eqn{n}. This function is very useful when you need to calculate the number of divisors for many numbers.
}
\value{
Returns an integer vector
}
\author{
Joseph Wood
}
\references{
\href{https://en.wikipedia.org/wiki/Divisor_function}{Divisor function}
}
\examples{
## Generate some random data
set.seed(8128)
mySamp <- sample(10^6, 5*10^5)

## Quickly generate number of coprime elements for many numbers
system.time(mySigmas <- numDivisorSieve(10^6))

## Now use result in algorithm
for (s in mySamp) {
    sSig <- mySigmas[s]
    ## Continue algorithm
}

}