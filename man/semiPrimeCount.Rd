\name{semiPrimeCount}
\encoding{UTF-8}
\alias{semiPrimeCount}
\title{
Count SemiPrimes With an Optional Partial-Sieve Cutoff
}
\description{
Count positive integers less than or equal to \eqn{n} with exactly two
prime factors, counted with multiplicity, excluding the first \eqn{a}
primes as factors. This is the correction term \eqn{P_2(n, a)} used in
Meissel and Lehmer prime counting.
}
\usage{
semiPrimeCount(n, a = 0, nThreads = NULL)
}
\arguments{
    \item{n}{A single nonnegative whole number, supplied as an integer or
    numeric value, at most \eqn{2^{53} - 1}{2^53 - 1}.}
    \item{a}{A single nonnegative whole number, supplied as an integer or
    numeric value, at most \eqn{2^{53} - 1}{2^53 - 1}. This is the number of
    initial primes excluded as factors, not a prime value. The default,
    \code{a = 0}, permits every prime. For example, \code{a = 1} excludes
    the prime 2, so only odd semiprimes are counted.}
    \item{nThreads}{Requested number of threads. The default is
    \code{NULL} for serial execution. Requests are capped by the available
    hardware and the useful work in the calculation. Small or simple cases
    may use fewer threads. This argument is validated even when the result
    is immediate.}
}
\details{
\subsection{Mathematical definition}{
Let \eqn{p_i} denote the \eqn{i}th prime. The function counts products
\eqn{p_i p_j \le n}{p[i] * p[j] <= n} with
\eqn{a < i \le j}{a < i <= j}. The ordering counts each integer once:
\eqn{2 \cdot 3}{2 * 3} and \eqn{3 \cdot 2}{3 * 2} do not give two
different semiprimes. Equal factors are allowed, so prime squares such as
4, 9, and 25 are included when their prime factor is permitted.

When \eqn{a = 0}, the excluded set is empty; no zeroth prime is needed.
The ten semiprimes at most 30 are
\code{4, 6, 9, 10, 14, 15, 21, 22, 25, 26}.
Excluding the first prime leaves \code{9, 15, 21, 25}, so
\code{semiPrimeCount(30, 1)} is 4.

Write \eqn{\pi(t)} for the number of primes at most \eqn{t}, and set
\eqn{b = \pi(\lfloor \sqrt n \rfloor)}{b = pi(floor(sqrt(n)))}. For a
fixed smaller factor \eqn{p_i}, the larger factor has an index from
\eqn{i} through \eqn{\pi(\lfloor n / p_i \rfloor)}{pi(floor(n / p[i]))}.
Consequently,
\deqn{P_2(n, a) = \sum_{i=a+1}^{b} \left[\pi(\lfloor n / p_i \rfloor) - i + 1\right].}{P2(n, a) = sum(i = a + 1, ..., b) [pi(floor(n / p[i])) - i + 1].}
A sum with no terms is zero. In particular, the answer is zero for
\eqn{n < 4} or \eqn{a \ge b}{a >= b}. Equivalently, the smallest
permitted prime square must be at most \eqn{n} for any semiprime to qualify.
}
\subsection{Relation to prime counting}{
For Meissel's cutoff
\eqn{a = \pi(\lfloor n^{1/3} \rfloor)}{a = pi(floor(n^(1/3)))},
the partial sieve leaves 1, the primes beyond the cutoff, and qualifying
semiprimes. Thus, for \eqn{n \ge 1}{n >= 1},
\deqn{\pi(n) = \Phi(n, a) + a - 1 - P_2(n, a).}{pi(n) = Phi(n, a) + a - 1 - P2(n, a).}
The function \code{\link{primePhi}} supplies \eqn{\Phi(n, a)}. Lehmer's
smaller cutoff also requires \code{\link{triPrimeCount}}; see
\code{\link{primeCount}} for the full formula.
}
}
\note{
The calculation prepares prime and prime-count tables with storage
proportional to \eqn{\sqrt n}. Larger prime-count queries in the sum use
the shared partial-sieve evaluator or a segmented sieve, selected according
to the input, cutoff, and useful thread budget. The segmented sieve uses
bounded buffers for each worker; the function does not allocate a dense
table through \eqn{n / p_{a+1}}{n / p[a + 1]}.

At the largest supported \eqn{n}, the square-root table has limit
94,906,265. Its dense prime-count entries alone require
\eqn{4(L + 1)}{4 * (L + 1)} bytes for limit \eqn{L}, about 362 MiB at
the maximum, excluding prime vectors, the shared Phi cache, temporary
storage, and allocation overhead. Workers within one call share the
prepared tables; separate calls prepare their own tables. A valid numeric
input may still require substantial time or memory.
}
\value{
A single whole number giving \eqn{P_2(n, a)}. The result has type
\code{integer} when it is at most \code{.Machine$integer.max}, and type
\code{double} otherwise. The count is exact throughout the supported range.
}
\author{
Joseph Wood
}
\seealso{
\code{\link{primePhi}}, \code{\link{triPrimeCount}},
\code{\link{primeCount}}, \code{\link{primeSieve}}
}
\examples{
## Count all semiprimes, then exclude 2 as a factor
semiPrimeCount(30)    ## 10
semiPrimeCount(30, 1) ## 4

## Equal prime factors are included
semiPrimeCount(4)     ## 1, from 2 * 2
semiPrimeCount(3)     ## 0
semiPrimeCount(0)     ## 0
semiPrimeCount(30, 100) ## 0

## A modest calculation with a requested thread budget
semiPrimeCount(10^6, nThreads = 2)
}
\keyword{prime}
\keyword{semiprime}
\keyword{number theory}
