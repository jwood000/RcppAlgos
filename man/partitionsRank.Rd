\name{partitionsRank}
\alias{partitionsRank}
\alias{compositionsRank}
\title{
Rank Partitions/Compositions
}
\description{
    \itemize{
        \item Generate the rank (lexicographically) of partitions/compositions. These functions are the complement to \code{partitions/compositionsSample}. See the examples below.
        \item GMP support allows for exploration of partitions/compositions of vectors with many elements.
    }
}
\usage{
partitionsRank(..., v, repetition = FALSE, freqs = NULL, target = NULL)

compositionsRank(..., v, repetition = FALSE, freqs = NULL,
                 target = NULL, weak = FALSE)
}
\arguments{
  \item{...}{vectors or matrices to be ranked.}
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. All atomic types are supported (See \code{\link{is.atomic}}).}
  \item{repetition}{Logical value indicating whether partitions/compositions should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all partitions of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{target}{Number to be partitioned. If \code{NULL}, the \code{max(v)} will be used.}
  \item{weak}{(Compositions only) Logical flag indicating whether to allow terms of the sequence to be zero.}
 }
\details{
These algorithms rely on efficiently ranking the \eqn{n^{th}}{nth} lexicographical partition.
}
\value{A vector of class \code{integer}, \code{numeric}, or \code{bigz} determined by the total number of partitions/compositions}
\references{
\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical order}
\href{https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation}{ranking/unranking}
}
\seealso{
\code{\link{partitionsSample}}, \code{\link{compositionsSample}}
}
\author{
Joseph Wood
}
\note{
    \code{v} must be supplied.
}
\examples{
mySamp = partitionsSample(30, 8, TRUE, n = 5, seed = 10, namedSample = TRUE)
myRank = partitionsRank(mySamp, v = 30, repetition = TRUE)
all.equal(as.integer(rownames(mySamp)), myRank)
}
\keyword{random}
\keyword{sample}
\keyword{rank}
\keyword{partitions}
\keyword{combinatorics}
