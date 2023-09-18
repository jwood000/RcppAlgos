\name{partitionsIter}
\alias{partitionsIter}
\alias{compositionsIter}
\title{
Partition/Composition Iterator
}
\description{
    \itemize{
        \item Returns an iterator for iterating over partitions/compositions of a numbers.
        \item Supports random access via the \code{[[} method.
        \item GMP support allows for exploration of cases where the number of partitions/compositions is large.
        \item Use the \code{next} methods to obtain results in lexicographical order.
    }
}
\usage{
partitionsIter(v, m = NULL, repetition = FALSE,
               freqs = NULL, target = NULL,
               nThreads = NULL, tolerance = NULL)

compositionsIter(v, m = NULL, repetition = FALSE, freqs = NULL,
                 target = NULL, weak = FALSE, nThreads = NULL,
                 tolerance = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. Only integer and numeric vectors are accepted.}
  \item{m}{Width of the partition. If \code{m = NULL}, the length will be determined by the partitioning case (\emph{e.g.} When we are generating distinct partitions of \eqn{n}, the width will be equal to the smallest \eqn{m} such that \code{sum(1:m) >= n}).}
  \item{repetition}{Logical value indicating whether partitions/compositions should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all partitions of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{target}{Number to be partitioned. If \code{NULL}, \code{max(v)} will be used.}
  \item{weak}{(Compositions only) Logical flag indicating whether to allow terms of the sequence to be zero.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
  \item{tolerance}{A numeric value greater than or equal to zero. This parameter is utilized when a constraint is applied on a numeric vector. The default value is 0 when it can be determined that whole values are being utilized, otherwise it is \code{sqrt(.Machine$double.eps)} which is approximately \eqn{1.5e-8}. N.B. If the input vector is of type integer, this parameter will be ignored and strict equality will be enforced.}
}

\value{
    \itemize{
        \item If \code{nextIter} is called, a vector is returned
        \item Otherwise, a matrix with \eqn{m} columns
    }
}
\details{
Once you initialize a new iterator, the following methods are available:
  \describe{
    \item{\code{nextIter}}{Retrieve the \bold{next} lexicographical result}
    \item{\code{nextNIter}}{Pass an integer \emph{n} to retrieve the \bold{next} \emph{n} lexicographical results}
    \item{\code{nextRemaining}}{Retrieve all remaining lexicographical results}
    \item{\code{currIter}}{Returns the current iteration}
    \item{\code{startOver}}{Resets the iterator}
    \item{\code{sourceVector}}{View the source vector}
    \item{\code{summary}}{Returns a list of summary information about the iterator}
    \item{\code{front}}{Retrieve the \bold{first} lexicographical result}
    \item{\code{back}}{Retrieve the \bold{last} lexicographical result}
    \item{\code{[[}}{Random access method. Pass a single value or a vector of valid indices. If a single value is passed, the internal index of the iterator will be updated, however if a vector is passed the internal state will not change. GMP support allows for flexible indexing.}
  }
}
\note{
    \itemize{
        \item{If \code{nThreads} is utilized, it will only take effect if the number of elements requested is greater than some threshold (determined internally). \emph{E.g}:
            \preformatted{serial   <- partitionsIter(1000, 10)
multi    <- partitionsIter(1000, 10, nThreads = 4)
fetch1e6 <- multi@nextNIter(1e6)  ## much faster than serial@nextNIter(1e6)
fetch1e3 <- multi@nextNIter(1e3)  ## only one thread used... same as serial@nextNIter(1e3)

library(microbenchmark)
microbenchmark(multi@nextNIter(1e6), serial@nextNIter(1e6))
microbenchmark(multi@nextNIter(1e3), serial@nextNIter(1e3))}
        }
        \item{\code{nThreads} will be ignored in the following cases (i.e. Generating the \eqn{n^{th}}{nth} partition in these cases are currently unavailable):
            \itemize{
                \item{With standard multisets. If zero is the only element with a non-trivial multiplicity, multithreading is possible.}
                \item{If the source vector is not isomorphic to \code{1:length(v)}}
            }
        }
        \item{The maximum number of partitions/compositions that can be generated at one time is \eqn{2^{31} - 1}{2^31 - 1}.}
    }
}
\seealso{
\code{\link{partitionsGeneral}}, \code{\link{compositionsGeneral}}
}
\references{
    \itemize{
        \item{\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical Order}}
        \item{\href{https://en.wikipedia.org/wiki/Subset_sum_problem}{Subset Sum Problem}}
        \item{\href{https://en.wikipedia.org/wiki/Partition_(number_theory)}{Partition (number theory)}}
        \item{\href{https://en.wikipedia.org/wiki/Composition_(combinatorics)}{Composition (combinatorics))}}
    }
}
\author{
Joseph Wood
}
\examples{
a = partitionsIter(0:10, repetition = TRUE)
a@nextIter()
a@nextNIter(3)
a@front()
a@nextRemaining()
a@summary()
a@back()
a[[5]]
a@summary()
a[[c(1, 17, 3)]]
a@summary()

## Multisets... no random access
b = partitionsIter(40, 5, freqs = rep(1:4, 10), target = 80)
b@nextIter()
b@nextNIter(10)
b@summary()
b@nextIter()
b@currIter()
}
\keyword{integer}
\keyword{partitions}
\keyword{combinatorics}
