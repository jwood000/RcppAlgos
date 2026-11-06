\name{primeCount}
\encoding{UTF-8}
\alias{primeCount}
\title{
Prime Counting Function \eqn{\pi(x)}{\pi(x)}
}
\description{
Count the prime numbers less than or equal to \eqn{n}, using automatic
selection or an explicitly selected Legendre, Meissel, or Lehmer method.
}
\usage{
primeCount(n, nThreads = NULL, method = "auto")
}
\arguments{
    \item{n}{A single positive whole number, at most
    \eqn{2^{53} - 1}{2^53 - 1}.}
    \item{nThreads}{A positive whole number requesting the number of threads,
    or \code{NULL} (the default) for serial execution. Requests are capped by
    the available hardware, and a method may use fewer threads. For
    \code{n < 100000}, this argument is ignored and execution is serial.}
    \item{method}{One of \code{"auto"}, \code{"legendre"},
    \code{"meissel"}, or \code{"lehmer"}. Matching is case-sensitive and
    accepts unique abbreviations, as in \code{\link{match.arg}}.
    \code{NULL} selects \code{"auto"}.}
}
\details{
\subsection{Method selection}{
The current \code{"auto"} policy uses small lookup values or a sieve for
\code{n < 100000}, and Legendre's formula otherwise. This preserves the
previous default behavior. The policy may change after further benchmarking;
\code{"auto"} does not currently guarantee the fastest method for each input.

An explicit method uses the named algorithm for every accepted \code{n},
including small inputs. The \code{method} argument follows \code{nThreads},
so existing positional calls such as \code{primeCount(n, 2)} retain their
meaning. All methods use serial execution when \code{nThreads = NULL}.
}
\subsection{Counting formulas}{
Let \eqn{p_a} be the \eqn{a}th prime, and let \eqn{\Phi(n, a)} count the
positive integers at most \eqn{n} that are not divisible by any of the first
\eqn{a} primes. Its recurrence is
\deqn{\Phi(n, a) = \Phi(n, a - 1) - \Phi(\lfloor n / p_a \rfloor, a - 1)}{Phi(n, a) = Phi(n, a - 1) - Phi(floor(n / p[a]), a - 1)}
with \eqn{\Phi(n, 0) = n}{Phi(n, 0) = n}.

For \eqn{k = 2} or \eqn{3}, let \eqn{P_k(n, a)} count integers at most
\eqn{n} with exactly \eqn{k} prime factors, none among the first
\eqn{a} primes. Factors are counted with multiplicity: \eqn{p^2} is included
in \eqn{P_2}, and \eqn{p^3} and \eqn{p^2 q} are included in \eqn{P_3} when
their factors exceed the cutoff. Nondecreasing factor order counts each
integer once. When \eqn{a = 0}, all primes are allowed.

The methods use different cutoffs and correction terms:
\itemize{
\item Legendre uses \eqn{a = \pi(\lfloor \sqrt n \rfloor)}{a = pi(floor(sqrt(n)))}:
\deqn{\pi(n) = \Phi(n, a) + a - 1.}{pi(n) = Phi(n, a) + a - 1.}
\item Meissel uses \eqn{a = \pi(\lfloor n^{1/3} \rfloor)}{a = pi(floor(n^(1/3)))}:
\deqn{\pi(n) = \Phi(n, a) + a - 1 - P_2(n, a).}{pi(n) = Phi(n, a) + a - 1 - P2(n, a).}
\item Lehmer uses \eqn{a = \pi(\lfloor n^{1/4} \rfloor)}{a = pi(floor(n^(1/4)))}:
\deqn{\pi(n) = \Phi(n, a) + a - 1 - P_2(n, a) - P_3(n, a).}{pi(n) = Phi(n, a) + a - 1 - P2(n, a) - P3(n, a).}
}

The Legendre implementation and shared \eqn{\Phi} calculation build on
Kim Walisch's work in \href{https://github.com/kimwalisch/primecount}{primecount},
including cached partial-sieve values, small-prime periodic tables, and
prime-count lookup shortcuts.
}
}
\note{
The counting methods prepare tables with storage proportional to
\eqn{\sqrt n}. At the maximum input \eqn{2^{53} - 1}{2^53 - 1}, the dense
prime-count table requires about 362 MiB and the prime-vector reservation
requests about 42 MiB, before caches and other overhead. These are calculated
storage estimates, not measured peak memory. Workers within one call share
the tables; concurrent calls prepare separate tables. The supported numeric
range does not guarantee that every input is practical for available memory.
}
\references{
\itemize{
\item Tomás Oliveira e Silva,
\href{https://sweet.ua.pt/tos/bib/5.4.pdf}{Computing \eqn{\pi(x)}: the combinatorial method},
Revista do DETUA, vol. 4, no. 6, March 2006, p. 761.
\item \href{https://en.wikipedia.org/wiki/Double-precision_floating-point_format}{53-bit significand precision}
}
}
\value{
A single whole number giving the number of primes less than or equal to
\eqn{n}. The result has type \code{integer} when it is at most
\code{.Machine$integer.max}, and type \code{double} otherwise. The count is
exact throughout the supported range.
}
\author{
Joseph Wood
}
\seealso{
\code{\link{primeSieve}}
}
\examples{
## Get the number of primes less than or equal to a billion
primeCount(10^9)

## Select a method explicitly
primeCount(10^6, method = "legendre")
primeCount(10^6, method = "meissel")
primeCount(10^6, method = "lehmer")

## Using nThreads
system.time(primeCount(10^10, nThreads = 2, method = "lehmer"))
}
\keyword{prime}
\keyword{pi}
\keyword{primecount}
\keyword{number theory}
