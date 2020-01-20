\name{isPrimeRcpp}
\alias{isPrimeRcpp}
\title{
Vectorized Primality Test
}
\description{
Implementation of the \href{https://en.wikipedia.org/wiki/Miller-Rabin_primality_test}{Miller-Rabin primality test}. Based on the "mp_prime_p" function from the "factorize.c" source file found in the gmp library: \url{http://gmplib.org}.
}
\usage{
isPrimeRcpp(v, namedVector = FALSE, nThreads = NULL)
}
\arguments{
    \item{v}{Vector of integers or numeric values.}
    \item{namedVector}{Logical flag. If \code{TRUE}, a named vector is returned. The default is \code{FALSE}.}
    \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
}
\note{
The maximum value for each element in \eqn{v} is \eqn{2^{53} - 1}{2^53 - 1}.
}
\details{
The Miller-Rabin primality test is a probabilistic algorithm that makes heavy use of \href{https://en.wikipedia.org/wiki/Modular_exponentiation}{modular exponentiation}. At the heart of modular exponentiation is the ability to accurately obtain the remainder of the product of two numbers \eqn{\pmod p}{mod p}.

With the gmp library, producing accurate calculations for problems like this is trivial because of the nature of the multiple precision data type. However, standard C++ does not afford this luxury and simply relying on a strict translation would have limited this algorithm to numbers less than \eqn{\sqrt 2^{63} - 1}{sqrt(2^63 - 1)} (N.B. We are taking advantage of the signed 64-bit fixed width integer from the stdint library in C++. If we were confined to base R, the limit would have been \eqn{\sqrt 2^{53} - 1}{sqrt(2^53 - 1)}). RcppAlgos::isPrimeRcpp gets around this limitation with a \href{https://en.wikipedia.org/wiki/Divide_and_conquer_algorithm}{divide and conquer} approach taking advantage of properties of arithmetic.

The problem we are trying to solve can be summarized as follows:

\deqn{(x_1 * x_2) \pmod p}{(x1 * x2) mod p}

Now, we rewrite \eqn{x_2}{x2} as \eqn{x_2 = y_1 + y_2 + \dots + y_n}{x2 = y1 + y2 + \dots + yn}, so that we obtain:

\deqn{(x_1 * y_1) \pmod p + (x_1 * y_2) \pmod p + \dots + (x_1 * y_n) \pmod p}{(x1 * y1) mod p + (x1 * y2) mod p + \dots + (x1 * yn) mod p}

Where each product \eqn{(x_1 * y_j)}{x1 * yj} for \eqn{j <= n} is smaller than the original \eqn{x_1 * x_2}{x1 * x2}. With this approach, we are now capable of handling much larger numbers. Many details have been omitted for clarity.

For a more in depth examination of this topic see \href{https://codereview.stackexchange.com/questions/186751/accurate-modular-arithmetic-with-double-precision}{Accurate Modular Arithmetic with Double Precision}.
}
\value{
Returns a named/unnamed logical vector. If an index is \code{TRUE}, the number at that index is prime, otherwise the number is composite.
}
\references{
\itemize{
\item{\href{http://www.math.uconn.edu/~kconrad/blurbs/ugradnumthy/millerrabin.pdf}{THE MILLER-RABIN TEST}\itemize{\item{Conrad, Keith. "THE MILLER-RABIN TEST." http://www.math.uconn.edu/~kconrad/blurbs/ugradnumthy/millerrabin.pdf.}}}
\item{\href{https://en.wikipedia.org/wiki/Double-precision_floating-point_format}{53-bit significand precision}}}
}
\seealso{
\code{\link{primeFactorize}}, \code{\link[gmp]{isprime}}, \code{\link[numbers]{isPrime}}
}
\examples{
## check the primality of a single number
isPrimeRcpp(100)

## check the primality of every number in a vector
isPrimeRcpp(1:100)

set.seed(42)
mySamp <- sample(10^13, 10)

## return named vector for easy identification
isPrimeRcpp(mySamp, namedVector = TRUE)

## Using nThreads
system.time(isPrimeRcpp(mySamp, nThreads = 2))
}