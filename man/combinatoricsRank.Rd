\name{comboRank}
\alias{comboRank}
\alias{permuteRank}
\title{
Rank Combinations and Permutations
}
\description{
    \itemize{
        \item Generate the rank (lexicographically) of combinations/permutations. These functions are the complement to \code{comboSample} and \code{permuteSample}. See the examples below.
        \item GMP support allows for exploration of combinations/permutations of vectors with many elements.
    }
}
\usage{
comboRank(..., v, repetition = FALSE, freqs = NULL)
permuteRank(..., v, repetition = FALSE, freqs = NULL)
}
\arguments{
  \item{...}{vectors or matrices to be ranked.}
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. All atomic types are supported (See \code{\link{is.atomic}}).}
  \item{repetition}{Logical value indicating whether combinations/permutations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all combinations/permutations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
 }
\details{
These algorithms rely on efficiently ranking the \eqn{n^{th}}{nth} lexicographical combination/permutation.
}
\value{A vector of class \code{integer}, \code{numeric}, or \code{bigz} determined by the total number of combinations/permutations}
\references{
\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical order}
\href{https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation}{ranking/unranking}
}
\seealso{
\code{\link{comboSample}}, \code{\link{permuteSample}}
}
\author{
Joseph Wood
}
\note{
    \code{v} must be supplied.
}
\examples{
mySamp = comboSample(30, 8, TRUE, n = 5, seed = 10, namedSample = TRUE)
myRank = comboRank(mySamp, v = 30, repetition = TRUE)
all.equal(as.integer(rownames(mySamp)), myRank)
}
\keyword{random}
\keyword{sample}
\keyword{rank}
\keyword{combinations}
\keyword{permutations}
\keyword{combinatorics}
