\name{primeFactorizationList}
\alias{primeFactorizationList}
\title{
Generate Prime Factorization for all Integers up to \eqn{n}
}
\description{
Efficiently generates the prime factorization of all integers up to \eqn{n} using Rcpp.
}
\usage{
primeFactorizationList(n = 100L)
}
\arguments{
  \item{n}{Positive integer or numeric value.}
}
\details{
This function is useful when many prime factorizations are needed. Instead of generating the prime factorization on the fly, one can reference the indices of the generated list. 
}
\value{
Returns a list of integer vectors where each index, \eqn{j}, represents the prime factorization of \eqn{j} (\emph{E.g.} \code{gmp::factorize(j)}).
}
\author{
Joseph Wood
}
\note{
The maximum allowed value is \eqn{2^{31} - 1}{2^31 - 1}.
}
\seealso{
\code{\link{divisorsList}}, \code{\link[gmp]{factorize}}, \code{\link[numbers]{primeFactors}}
}
\examples{
## Generate some random data
set.seed(28)
mySamp <- sample(10^6, 5*10^5)

## Quickly generate prime factorizations up 
## to 10^6 (max element from mySamp)
system.time(allPFacs <- primeFactorizationList(10^6))

## Use generated prime factorization for further
## analysis by accessing the index of allPFacs
for (s in mySamp) {
    pFac <- allPFacs[[s]]
    ## Continue algorithm
}

## Generating the prime factorization on the fly is much 
## slower despite 1. gmp::factorize being one of the fastest
## stand alone factoring algorithms and 2. there are fewer
## factorizations being performed (allPFacs above contains
## 10^6 factorizations and the below only generates the
## factorizations of the numbers in mySamp (5*10^5 elements)).

\dontrun{
system.time(mySampPFacs <- lapply(mySamp, gmp::factorize))}
}
