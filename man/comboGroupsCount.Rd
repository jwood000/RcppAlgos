\name{comboGroupsCount}
\alias{comboGroupsCount}
\title{
Number of Partitions of a Vector into Groups
}
\description{
Calculate the number of partitions of a vector into groups. See the related integer sequences A025035-A025042 at \href{https://oeis.org}{OEIS} (E.g. \href{https://oeis.org/A025036}{A025036} for Number of partitions of { 1, 2, ..., 4n } into sets of size 4.)
}
\usage{
comboGroupsCount(v, numGroups = NULL, grpSizes = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. All atomic types are supported (See \code{\link{is.atomic}}).}
  \item{numGroups}{An Integer. The number of groups that the vector will be partitioned into. The default is \code{NULL}. If provided and \code{grpSize} is \code{NULL}, it must divide the length of v (if v is a vector) or v (if v is a scalar).}
  \item{grpSizes}{A vector of whole numbers representing the size of each group. The default is \code{NULL}. If provided, the sum of the elements must total the length of v (if v is a vector) or v (if v is a scalar).}
}
\value{
A numerical value representing the total number of partitions of groups.
}
\references{
\href{https://oeis.org/A025036}{OEIS Integer Sequence A025036}
}
\author{
Joseph Wood
}
\note{
When the number of results exceeds \eqn{2^{53} - 1}{2^53 - 1}, a number of class \code{bigz} is returned.
}
\examples{
comboGroupsCount(16, 4)
comboGroupsCount(16, grpSizes = c(1:4, 6))
comboGroupsCount(28, grpSizes = rep(2:5, each = 2))
}
\keyword{partitions}
\keyword{comboGroups}
\keyword{combinatorics}
