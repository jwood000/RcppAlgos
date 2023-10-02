\name{partitionsSample}
\alias{partitionsSample}
\alias{compositionsSample}
\title{
Sample Partitions/Compositions
}
\description{
    \itemize{
        \item Generate a specific (lexicographically) or random sample of partitions/compositions of a number.
        \item Produce results in parallel using the \code{Parallel} or \code{nThreads} arguments.
        \item GMP support allows for exploration of cases where the number of partitions/compositions is large.
    }
}
\usage{
partitionsSample(v, m = NULL, repetition = FALSE, freqs = NULL,
                 target = NULL, n = NULL, sampleVec = NULL,
                 seed = NULL, nThreads = NULL, namedSample = FALSE)

compositionsSample(v, m = NULL, repetition = FALSE, freqs = NULL,
                   target = NULL, weak = FALSE, n = NULL, sampleVec = NULL,
                   seed = NULL, nThreads = NULL, namedSample = FALSE)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. Only integer and numeric vectors are accepted.}
  \item{m}{Width of the partition. If \code{m = NULL}, the length will be determined by the partitioning case (\emph{e.g.} When we are generating distinct partitions of \eqn{n}, the width will be equal to the smallest \eqn{m} such that \code{sum(1:m) >= n}).}
  \item{repetition}{Logical value indicating whether partitions/compositions should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all partitions of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{target}{Number to be partitioned. If \code{NULL}, \code{max(v)} will be used.}
  \item{weak}{(Compositions only) Logical flag indicating whether to allow terms of the sequence to be zero.}
  \item{n}{Number of partitions/compositions to return. The default is \code{NULL}.}
  \item{sampleVec}{A vector of numbers representing the lexicographical partitions/compositions to return. Accepts vectors of class \code{bigz} as well as vectors of characters}
  \item{seed}{Random seed initialization. The default is \code{NULL}. N.B. If the gmp library is needed, this parameter must be set in order to have reproducible results (\emph{E.g} \code{set.seed()} has no effect in these cases).}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}.}
  \item{namedSample}{Logical flag. If \code{TRUE}, \code{rownames} corresponding to the lexicographical partition, will be added to the returned matrix. The default is \code{FALSE}.}
 }
\details{
These algorithms rely on efficiently generating the \eqn{n^{th}}{nth} lexicographical partition. This is the process of \href{https://www.rosettacode.org/wiki/Permutations/Rank_of_a_permutation}{unranking}.
}
\value{
A matrix is returned with each row containing a vector of length \eqn{m}.
}
\references{
    \itemize{
        \item{\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical order}}
        \item{\href{https://en.wikipedia.org/wiki/Partition_(number_theory)}{Partition (Number Theory)}}
    }
}
\author{
Joseph Wood
}
\note{
    \itemize{
        \item{\code{partitionsSample} is not available for the following cases:
            \itemize{
                \item With standard multisets. If zero is the only element with a non-trivial multiplicity, sampling is allowed (\emph{e.g.} \code{partitionsSample(0:100, freqs = c(100, rep(1, 100)), n = 2)})
                \item If the source vector is not isomorphic to \code{1:length(v)} (\emph{e.g.} \code{v = c(1, 4, 6, 7, 8)}).
            }
        }
        \item \code{n} and \code{sampleVec} cannot both be \code{NULL}.
    }
}

\examples{
partitionsSample(100, 10, n = 5)
partitionsSample(100, 10, seed = 42, n = 5, target = 200)

## retrieve specific results (lexicographically)
partitionsCount(100, 10, TRUE, target = 500)
## [1] 175591757896
partitionsSample(100, 10, TRUE, target = 500,
                 sampleVec = c(1, 1000, 175591757896))
}
\keyword{random}
\keyword{sample}
\keyword{rank}
\keyword{integer}
\keyword{unrank}
\keyword{partitions}
\keyword{combinatorics}
