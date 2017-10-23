\name{divisorsList}
\alias{divisorsList}
\title{
Generate Complete Factorization of all Integers up to \eqn{n}
}
\description{
Rcpp sieve implementation that quickly generates the complete factorization of all integers up to \eqn{n}.
}
\usage{
divisorsList(n = 100L)
}
\arguments{
  \item{n}{Positive integer or numeric value.}
}
\details{
This function is useful when many complete factorizations are needed. Instead of generating the complete factorization on the fly, one can reference the indices of the generated list. 
}
\value{
Returns a list of integer vectors where each index, \eqn{m}, represents the complete factorization of \eqn{m}.

\emph{E.g.} \code{(1:m)[m\%\%(1:m)==0]}
}
\author{
Joseph Wood
}
\seealso{
\code{\link{primeFactorizationList}}, \code{\link[numbers]{divisors}}
}
\examples{
## Generate some random data
set.seed(33550336)
mySamp <- sample(10^5, 5*10^4)

## Quickly generate complete factorizations up 
## to 10^5 (max element from mySamp)
system.time(allFacs <- divisorsList(10^5))

## Use generated complete factorization for further
## analysis by accessing the index of allFacs
for (s in mySamp) {
    myFac <- allFacs[[s]]
    ## Continue algorithm
}

## Generating the complete factorization on the fly is much slower.
## See below for a couple of examples that are not as efficient.

\dontrun{
## Brute Force... very slow
system.time(mySampFacs <- lapply(mySamp, function(x) (1:x)[x\%\%(1:x)==0L]))

## Using numbers::divisors... better than brute force, but still
## much slower than calling divisorList
system.time(mySampFacs <- lapply(mySamp, function(x) numbers::divisors))
}
}