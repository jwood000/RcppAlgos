\name{primeSieve}
\alias{primeSieve}
\title{
Generate Prime Integers
}
\description{
Rcpp implementation of the sieve of Eratosthenes. Generates all prime numbers between \code{bound1} and \code{bound2} (if supplied) or all primes up to \code{bound1}.
}
\usage{
primeSieve(bound1 = 100L, bound2 = NULL)
}
\arguments{
\item{bound1}{Positive integer or numeric value.}

\item{bound2}{Positive integer or numeric value.}
}
\note{
It does not matter which bound is larger as the resulting primes will be between \code{min(bound1, bound2)} and \code{max(bound1, bound2)} if \code{bound2} is provided.
}
\details{
Optimized implementation of the sieve of Eratosthenes. In addition to being very fast, this algorithm is particularly memory efficient when the sieving interval isn't large (\emph{E.g.} check \code{object.size} when \code{abs(bound1 - bound2) < 10^5}). This is accomplished by doing most of the work with logical vectors and only allocating enough memory for the return vector based off of the size of the interval.
}
\value{
Returns an integer vector
}
\seealso{
\code{\link[numbers]{Primes}}
}
\references{
\href{https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes}{Sieve of Eratosthenes}
}
\examples{
## Primes up to a thousand
primeSieve(1000)

## Primes between 42 and 1729
primeSieve(42, 1729)

## Equivalent to 
primeSieve(1729, 42)

## Primes up to one hundred million in no time
system.time(primeSieve(10^8))

## Quickly generate large primes over small interval
system.time(myPs <- primeSieve(1000001000, 10^9))
## Object created is small
object.size(myPs)

## Compared to generating all primes up to the 
## upper bound and subsetting (see below)

\dontrun{
## Takes a decent amount of time
system.time(p <- primeSieve(1000001000))
system.time(p[p >= 10^9])
## Creates large object
object.size(p)}
}