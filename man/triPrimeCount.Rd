\name{triPrimeCount}
\encoding{UTF-8}
\alias{triPrimeCount}
\title{
Count TriPrimes With an Optional Partial-Sieve Cutoff
}
\description{
Count positive integers less than or equal to \eqn{n} with exactly three
prime factors, counted with multiplicity, excluding the first \eqn{a}
primes as factors. This is the correction term \eqn{P_3(n, a)} used in
Lehmer prime counting. Its table requirements impose an additional
input-dependent limit described below.
}
\usage{
triPrimeCount(n, a = 0, nThreads = NULL)
}
\arguments{
    \item{n}{A single nonnegative whole number, supplied as an integer or
    numeric value, at most \eqn{2^{53} - 1}{2^53 - 1}. The additional
    table limit in \sQuote{Details} may restrict the accepted value for a
    given \code{a}. With the default \code{a = 0}, the largest accepted
    \code{n} is 379,625,063.}
    \item{a}{A single nonnegative whole number, supplied as an integer or
    numeric value, at most \eqn{2^{53} - 1}{2^53 - 1}. This is the number of
    initial primes excluded as factors, not a prime value. The default,
    \code{a = 0}, permits every prime. For example, \code{a = 1} excludes
    the prime 2, so only odd triprimes are counted. Larger cutoffs can
    permit larger \code{n} within the table limit.}
    \item{nThreads}{A positive whole number at most
    \code{.Machine$integer.max}, or \code{NULL} (the default) for serial
    execution. Requests are capped by the available hardware and the useful
    work in the calculation. Small or simple cases may use fewer threads.
    This argument is validated even when the result is immediate.}
}
\details{
\subsection{Mathematical definition}{
Let \eqn{p_i} denote the \eqn{i}th prime. The function counts products
\eqn{p_i p_j p_k \le n}{p[i] * p[j] * p[k] <= n} with
\eqn{a < i \le j \le k}{a < i <= j <= k}. Nondecreasing factor order
counts each integer once. Factors need not be distinct: prime cubes
\eqn{p^3} and products \eqn{p^2 q} are included, as are products of
three different primes.

When \eqn{a = 0}, the excluded set is empty; no zeroth prime is needed.
The seven triprimes at most 30 are \code{8, 12, 18, 20, 27, 28, 30}.
For example, 8 qualifies as \eqn{2^3}, while 16 does not qualify because
\eqn{2^4} has four prime factors counted with multiplicity. Excluding 2
leaves only 27 among these seven values.

Write \eqn{\pi(t)} for the number of primes at most \eqn{t}, and define
\deqn{c = \pi(\lfloor n^{1/3} \rfloor), \qquad b_i = \pi(\lfloor \sqrt{n / p_i} \rfloor).}{c = pi(floor(n^(1/3))); b[i] = pi(floor(sqrt(n / p[i]))).}
The smallest factor has index at most \eqn{c}; for that factor, the middle
factor has index at most \eqn{b_i}. Once these two factors are fixed, the
largest factor has index from \eqn{j} through
\eqn{\pi(\lfloor n / (p_i p_j) \rfloor)}{pi(floor(n / (p[i] * p[j])))}.
Hence,
\deqn{P_3(n, a) = \sum_{i=a+1}^{c} \sum_{j=i}^{b_i} \left[\pi(\lfloor n / (p_i p_j) \rfloor) - j + 1\right].}{P3(n, a) = sum(i = a + 1, ..., c) sum(j = i, ..., b[i]) [pi(floor(n / (p[i] * p[j]))) - j + 1].}
A sum with no terms is zero. In particular, the answer is zero for
\eqn{n < 8} or \eqn{a \ge c}{a >= c}: the smallest permitted prime
cube then exceeds \eqn{n}.
}
\subsection{Relation to prime counting}{
At Lehmer's cutoff
\eqn{a = \pi(\lfloor n^{1/4} \rfloor)}{a = pi(floor(n^(1/4)))},
the surviving composites have two or three permitted prime factors.
For \eqn{n \ge 1}{n >= 1},
\deqn{\pi(n) = \Phi(n, a) + a - 1 - P_2(n, a) - P_3(n, a).}{pi(n) = Phi(n, a) + a - 1 - P2(n, a) - P3(n, a).}
The other terms are available through \code{\link{primePhi}} and
\code{\link{semiPrimeCount}}. See \code{\link{primeCount}} for the
complete counting methods.
}
\subsection{Table limit and memory use}{
The implementation first prepares tables through
\eqn{\lfloor n^{1/3} \rfloor}{floor(n^(1/3))} to determine whether any
qualifying triples exist. If they do, let \eqn{p_{a+1}} be the smallest
permitted prime. Direct prime-count lookups require a table through
\deqn{q_{\max} = \lfloor n / p_{a+1}^2 \rfloor.}{qmax = floor(n / p[a + 1]^2).}
Every queried quotient is at most this value. The largest supported table
limit is 94,906,265; a nonempty calculation requiring a larger table raises
an error before that larger table is allocated. An empty count does not
require this second table.

With \code{a = 0}, the smallest permitted prime is 2, so the condition is
\eqn{\lfloor n / 4 \rfloor \le 94906265}{floor(n / 4) <= 94906265}.
This explains the default upper limit \eqn{379625063}. Excluding more
primes can reduce the required table size and allow larger inputs, up to
the common numeric limit \eqn{2^{53} - 1}{2^53 - 1} when the table
condition is met.

Table storage is proportional to
\eqn{\max(n^{1/3}, n / p_{a+1}^2)}{max(n^(1/3), n / p[a + 1]^2)} for
nonempty calculations. The dense prime-count entries alone use
\eqn{4(L + 1)}{4 * (L + 1)} bytes for limit \eqn{L}, about 362 MiB at
the maximum. This excludes prime vectors, temporary tables, worker storage,
allocation overhead, and any shared Phi cache retained from other calls;
it is not a peak-memory bound. Workers within
one call share the prepared lookup tables. Inputs below the enforced table
limit can still be expensive in time or memory.
}
}
\value{
A single whole number giving \eqn{P_3(n, a)}. The result has type
\code{integer} when it is at most \code{.Machine$integer.max}, and type
\code{double} otherwise. The count is exact for accepted inputs.
}
\author{
Joseph Wood
}
\seealso{
\code{\link{primePhi}}, \code{\link{semiPrimeCount}},
\code{\link{primeCount}}, \code{\link{primeSieve}}
}
\examples{
## Count all triprimes, then exclude 2 as a factor
triPrimeCount(30)    ## 7
triPrimeCount(30, 1) ## 1, from 3 * 3 * 3

## Repeated prime factors are included
triPrimeCount(8)      ## 1, from 2 * 2 * 2
triPrimeCount(7)      ## 0
triPrimeCount(0)      ## 0
triPrimeCount(30, 100) ## 0

## Keep exploration modest: the default table extends through floor(n / 4)
triPrimeCount(10^5, nThreads = 2)
}
\keyword{prime}
\keyword{triprime}
\keyword{number theory}
