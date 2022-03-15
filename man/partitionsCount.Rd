\name{partitionsCount}
\alias{partitionsCount}
\title{
Number of Partitions
}
\description{
Calculate the number of partitions of a vector chosen \eqn{m} at a time with or without replacement. Additionally, these functions can calculate the number of partitions of multisets.
}
\usage{
partitionsCount(v, m = NULL, repetition = FALSE, freqs = NULL, target = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. Only integer and numeric vectors are accepted.}
  \item{m}{Width of the partition. If \code{m = NULL}, the length will be determined by the partitioning case (\emph{e.g.} When we are generating distinct partitions of \eqn{n}, the width will be equal to the smallest \eqn{m} such that \code{sum(1:m) >= n}).}
  \item{repetition}{Logical value indicating whether partitions should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all partitions of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{target}{Number to be partitioned.}
}
\value{
A numerical value representing the total number of partitions.
}
\note{
When the number of results exceeds \eqn{2^{53} - 1}{2^53 - 1}, a number of class \code{bigz} is returned.
}
\seealso{
\code{\link{partitionsGeneral}}
}
\examples{
## Same interface as partitionsGeneral
partitionsCount(25, 5)
partitionsCount(15, 7, TRUE)
partitionsCount(25, 5, freqs = rep(2, 25))

## Return object of class 'bigz'
partitionsCount(2500, 15, TRUE)
}
