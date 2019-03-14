\name{primeCount}
\alias{primeCount}
\title{
Prime Counting Function \eqn{\pi(x)}{\pi(x)}
}
\description{
\href{https://en.wikipedia.org/wiki/Prime-counting_function}{Prime counting function} for counting the prime numbers less than an integer, \eqn{n}, using Legendre's formula. It is based on the the algorithm developed by Kim Walisch found here: \href{https://github.com/kimwalisch/primecount}{kimwalisch/primecount}.
}
\usage{
primeCount(n, nThreads = NULL)
}
\arguments{
    \item{n}{Positive number}
    \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
}
\details{
\href{http://mathworld.wolfram.com/LegendresFormula.html}{Legendre's Formula} for counting the number of primes less than \eqn{n} makes use of the \href{https://en.wikipedia.org/wiki/Inclusion-exclusion_principle}{inclusion-exclusion principle} to avoid explicitly counting every prime up to \eqn{n}. It is given by:
\deqn{\pi(x) = \pi(\sqrt x) + \Phi(x, \sqrt x) - 1}
Where \eqn{\Phi(x, a)} is the number of positive integers less than or equal to \eqn{x} that are relatively prime to the first \eqn{a} primes (i.e. not divisible by any of the first \eqn{a} primes). It is given by the recurrence relation (\eqn{p_a}{pa} is the \eqn{ath} prime (e.g. \eqn{p_4 = 7}{p4 = 7})):
\deqn{\Phi(x, a) = \Phi(x, a - 1) + \Phi(x / p_a, a - 1)}{\Phi(x, a) = \Phi(x, a - 1) + \Phi(x / pa, a - 1)}
This algorithm implements five modifications developed by Kim Walisch for calculating \eqn{\Phi(x, a)} efficiently.
\enumerate{
\item{Cache results of \eqn{\Phi(x, a)}}
\item{Calculate \eqn{\Phi(x, a)} using \eqn{\Phi(x, a) = (x / pp) * \phi(pp) + \Phi(x mod pp, a)}  if \eqn{a <= 6}
\itemize{
\item{\eqn{pp = 2 * 3 * ... * } \code{prime[a]}}
\item{\eqn{\phi(pp) = (2 - 1) * (3 - 1) * ... * } \eqn{(}\code{prime[a]} \eqn{- 1)} (i.e. Euler's totient function)}
}
}
\item{Calculate \eqn{\Phi(x, a)} using \eqn{\pi(x)} lookup table}
\item{Calculate all \eqn{\Phi(x, a) = 1} upfront}
\item{Stop recursion at \eqn{6} if \eqn{\sqrt x >= 13} or \eqn{\pi(\sqrt x)} instead of \eqn{1}}
}
}
\note{
The maximum value of \eqn{n} is \eqn{2^{53} - 1}{2^53 - 1}
}
\references{
\itemize{
\item{\href{http://sweet.ua.pt/tos/bib/5.4.pdf}{Computing \eqn{\pi(x)}: the combinatorial method}\itemize{\item{Tomás Oliveira e Silva, Computing pi(x): the combinatorial method, Revista do DETUA, vol. 4, no. 6, March 2006, p. 761. http://sweet.ua.pt/tos/bib/5.4.pdf}}}
\item{\href{https://en.wikipedia.org/wiki/Double-precision_floating-point_format}{53-bit significand precision}}
}
}
\value{
Whole number representing the number of prime numbers less than or equal to \eqn{n}.
}
\author{
Joseph Wood
}
\seealso{
\code{\link{primeSieve}}
}
\examples{
## Get the number of primes less than a billion
primeCount(10^9)

## Using nThreads
system.time(primeCount(10^10, nThreads = 2))
}