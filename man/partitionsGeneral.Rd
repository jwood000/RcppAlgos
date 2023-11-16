\name{partitionsGeneral}
\alias{partitionsGeneral}
\alias{partitionsGeneral.default}
\alias{partitionsGeneral.table}
\alias{compositionsGeneral}
\alias{compositionsGeneral.default}
\alias{compositionsGeneral.table}
\title{
Generate Partitions/Compositions
}
\description{
The algorithms in \code{RcppAlgos} go beyond the traditional integer partition algorithms and can tackle a wide variety of cases.

    \itemize{
        \item Efficient algorithms for partitioning numbers under various constraints:
        \itemize{
            \item Standard (with repetition)
            \item Distinct
            \item Restricted
            \item Where each part has a specific multiplicity (i.e. when using \code{freqs} for multisets).
            \item Arbitrary target and source vector (\emph{e.g.} \code{partitionsGeneral(sample(1000, 20), 10, TRUE, target = 5000)})
        }
        \item Produce results in parallel using the \code{nThreads} arguments.
        \item Alternatively, the arguments \code{lower} and \code{upper} make it possible to generate partitions/compositions in chunks allowing for parallelization via the parallel package.
        \item GMP support allows for exploration of cases where the number of partitions/compositions is large.
        \item The output is in lexicographical order.
    }
}
\usage{
partitionsGeneral(v, m = NULL, ...)
compositionsGeneral(v, m = NULL, ...)

\method{partitionsGeneral}{default}(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL,
    lower = NULL, upper = NULL, nThreads = NULL, tolerance = NULL, ...
)
\method{compositionsGeneral}{default}(
    v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL, weak = FALSE,
    lower = NULL, upper = NULL, nThreads = NULL, tolerance = NULL, ...
)

\method{partitionsGeneral}{table}(
    v, m = NULL, target = NULL, lower = NULL,
    upper = NULL, nThreads = NULL, tolerance = NULL, ...
)
\method{compositionsGeneral}{table}(
    v, m = NULL, target = NULL, weak = FALSE, lower = NULL,
    upper = NULL, nThreads = NULL, tolerance = NULL, ...
)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. Only integer and numeric vectors are accepted.}
  \item{m}{Width of the partition. If \code{m = NULL}, the length will be determined by the partitioning case (\emph{e.g.} When we are generating distinct partitions of \eqn{n}, the width will be equal to the smallest \eqn{m} such that \code{sum(1:m) >= n}).}
  \item{\dots}{Further arguments passed to methods.}
  \item{repetition}{Logical value indicating whether partitions/compositions should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all partitions of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{lower}{The lower bound. Partitions/compositions are generated lexicographically, thus utilizing this argument will determine which specific partition to start generating from (\emph{e.g.} \code{partitionsGeneral(15, 3, lower = 6)} is equivalent to \code{partitionsGeneral(15, 3)[6:partitionsCount(15, 3), ]}). This argument along with \code{upper} is very useful for generating partitions/compositions in chunks allowing for easy parallelization.}
  \item{upper}{The upper bound. Similar to \code{lower}, however this parameter allows the user to \emph{stop} generation at a specific partition (\emph{e.g.} \code{partitionsGeneral(15, 3, upper = 5)} is equivalent to \code{partitionsGeneral(15, 3)[1:5, ]})}
  \item{target}{Number to be partitioned. If \code{NULL}, \code{max(v)} will be used.}
  \item{weak}{(Compositions only) Logical flag indicating whether to allow terms of the sequence to be zero.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
  \item{tolerance}{A numeric value greater than or equal to zero. This parameter is utilized when a constraint is applied on a numeric vector. The default value is 0 when it can be determined that whole values are being utilized, otherwise it is \code{sqrt(.Machine$double.eps)} which is approximately \eqn{1.5e-8}. N.B. If the input vector is of type integer, this parameter will be ignored and strict equality will be enforced.}
}
\value{
A matrix is returned with each row containing a vector of length \eqn{m}.
}
\note{
    \itemize{
        \item{\code{nThreads} will be ignored in the following cases (i.e. Generating the \eqn{n^{th}}{nth} partition in these cases are currently unavailable):
            \itemize{
                \item With standard multisets. If zero is the only element with a non-trivial multiplicity, multithreading is possible (\emph{e.g.} \code{partitionsGeneral(0:100, freqs = c(100, rep(1, 100)), nThreads = 4)}).
                \item If the source vector is not isomorphic to \code{1:length(v)} (\emph{e.g.} \code{v = c(1, 4, 6, 7, 8)}).
            }
        }
        \item The maximum number of partitions/compositions that can be generated at one time is \eqn{2^{31} - 1}{2^31 - 1}. Utilizing \code{lower} and \code{upper} makes it possible to generate additional partitions/compositions.
    }
}
\author{
Joseph Wood
}
\references{
    \itemize{
        \item{\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical Order}}
        \item{\href{https://en.wikipedia.org/wiki/Subset_sum_problem}{Subset Sum Problem}}
        \item{\href{https://en.wikipedia.org/wiki/Partition_(number_theory)}{Partition (number theory)}}
        \item{\href{https://en.wikipedia.org/wiki/Composition_(combinatorics)}{Composition (combinatorics))}}
    }
}
\examples{
partitionsGeneral(1)
partitionsGeneral(-1:0, 1)
partitionsGeneral(-1:0, 1, target = -1)
partitionsGeneral(20, 5)
partitionsGeneral(20, 5, repetition = TRUE)
partitionsGeneral(20, 5, freqs = rep(1:4, 5))
partitionsGeneral(20, 5, TRUE, target = 80)
partitionsGeneral(0:10, repetition = TRUE)
partitionsGeneral(seq(2L, 500L, 23L), 5, target = 1804)

compositionsGeneral(0:10, 5, repetition = TRUE)

set.seed(111)
partitionsGeneral(sample(1000, 20), 5, TRUE, target = 2500)

system.time(one_thread  <- partitionsGeneral(80, 10, TRUE))
system.time(two_threads <- partitionsGeneral(80, 10, TRUE, nThreads = 2))
identical(one_thread, two_threads)
}
\keyword{integer}
\keyword{partitions}
\keyword{combinatorics}
