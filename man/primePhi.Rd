\name{primePhi}
\encoding{UTF-8}
\alias{primePhi}
\title{
Count Integers Remaining After a Partial Prime Sieve
}
\description{
Count the positive integers less than or equal to \eqn{n} that are not
divisible by any of the first \eqn{a} primes. This is the partial-sieve
function \eqn{\Phi(n, a)} used in combinatorial prime counting, not
Euler's totient function.
}
\usage{
primePhi(n, a, nThreads = NULL)
}
\arguments{
    \item{n}{A single nonnegative whole number, supplied as an integer or
    numeric value, at most \eqn{2^{53} - 1}{2^53 - 1}.}
    \item{a}{A single nonnegative whole number, supplied as an integer or
    numeric value, at most \eqn{2^{53} - 1}{2^53 - 1}. This is the number of
    initial primes to exclude, not a prime value or an upper bound on the
    excluded primes. For example, \code{a = 3} excludes multiples of 2, 3,
    and 5. With \code{a = 0}, no primes are excluded.}
    \item{nThreads}{Requested number of threads. The default is
    \code{NULL} for serial execution. Requests are capped by the available
    hardware and the useful work in the calculation. Small or simple cases
    may use fewer threads. This argument is validated even when the result
    is immediate.}
}
\details{
\subsection{Mathematical definition}{
Write \eqn{p_1 = 2, p_2 = 3, p_3 = 5, \ldots}{p[1] = 2, p[2] = 3, p[3] = 5, ...}
for the increasing sequence of primes. Then \eqn{\Phi(n, a)} counts
integers \eqn{m} with \eqn{1 \le m \le n}{1 <= m <= n} for which none of
\eqn{p_1, \ldots, p_a}{p[1], ..., p[a]} divides \eqn{m}. When \eqn{a = 0},
the excluded set is empty. In particular,
\deqn{\Phi(n, 0) = n, \qquad \Phi(0, a) = 0.}{Phi(n, 0) = n; Phi(0, a) = 0.}
For \eqn{a \ge 1}{a >= 1}, removing multiples of the next prime gives
\deqn{\Phi(n, a) = \Phi(n, a - 1) - \Phi(\lfloor n / p_a \rfloor, a - 1).}{Phi(n, a) = Phi(n, a - 1) - Phi(floor(n / p[a]), a - 1).}
The subtracted term counts the multiples of \eqn{p_a} that survived
removal of the preceding primes.

The integer 1 survives every partial sieve. Thus, when \eqn{n > 0}, the
answer is at least 1, and it is exactly 1 when
\eqn{a \ge \pi(n)}{a >= pi(n)}, where \eqn{\pi(n)} is the number of primes
at most \eqn{n}. Composite numbers can also survive: for example,
\code{primePhi(30, 2)} includes 25 because neither 2 nor 3 divides it.
This gives the ten survivors
\code{1, 5, 7, 11, 13, 17, 19, 23, 25, 29}.

Euler's totient \eqn{\phi(n)} instead counts integers at most \eqn{n}
that are coprime to \eqn{n} itself. For that function, see
\code{\link{eulerPhiSieve}}.
}
\subsection{Relation to prime counting}{
If \eqn{a = \pi(\lfloor \sqrt n \rfloor)}{a = pi(floor(sqrt(n)))} and
\eqn{n \ge 1}{n >= 1}, every surviving integer other than 1 is a prime
larger than the excluded primes. Legendre's formula is therefore
\deqn{\pi(n) = \Phi(n, a) + a - 1.}{pi(n) = Phi(n, a) + a - 1.}
Smaller choices of \eqn{a} leave composite survivors. Meissel and Lehmer
prime counting account for these using the correction terms exposed by
\code{\link{semiPrimeCount}} and \code{\link{triPrimeCount}}. See
\code{\link{primeCount}} for the corresponding formulas.
}
}
\note{
The general calculation prepares prime and prime-count tables with storage
proportional to \eqn{\sqrt n}. Cases with \eqn{a \le 6}{a <= 6} use small
periodic tables instead of a full square-root table. Large values of \eqn{a}
do not require preparing all of the first \eqn{a} primes merely to establish
that the answer is 1.

The largest supported table limit is 94,906,265. Its dense prime-count table
alone requires \eqn{4(L + 1)}{4 * (L + 1)} bytes for limit \eqn{L}, about
362 MiB at that maximum. Prime vectors, the shared Phi cache, temporary
storage, and allocation overhead require additional memory. Workers within
one call share the prepared tables; separate calls prepare their own tables.
The numeric input range does not guarantee that every calculation is
practical for the available resources.
}
\value{
A single whole number giving \eqn{\Phi(n, a)}. The result has type
\code{integer} when it is at most \code{.Machine$integer.max}, and type
\code{double} otherwise. The count is exact throughout the supported range.
}
\author{
Joseph Wood
}
\seealso{
\code{\link{primeCount}}, \code{\link{semiPrimeCount}},
\code{\link{triPrimeCount}}, \code{\link{eulerPhiSieve}},
\code{\link{primeSieve}}
}
\examples{
## Remove multiples of the first two primes, 2 and 3
primePhi(30, 2)  ## 10

## No excluded primes; an empty interval; all primes through n excluded
primePhi(30, 0)   ## 30
primePhi(0, 2)    ## 0
primePhi(30, 100) ## 1

## Recover the prime count using Legendre's formula
n <- 1000
a <- primeCount(floor(sqrt(n)))
primePhi(n, a) + a - 1 == primeCount(n)
}
\keyword{prime}
\keyword{phi}
\keyword{sieve}
\keyword{number theory}
