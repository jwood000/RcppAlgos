\name{comboGroupsCount}
\alias{comboGroupsCount}
\title{
Number of Partitions of a Vector into Groups of Equal Size
}
\description{
Calculate the number of partitions of a vector into groups of equal size. See the related integer sequences A025035-A025042 at \href{https://oeis.org}{OEIS} (E.g. \href{https://oeis.org/A025036}{A025036} for Number of partitions of { 1, 2, ..., 4n } into sets of size 4.)
}
\usage{
comboGroupsCount(v, numGroups)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{Source vector. If \code{v} is an integer (including nonpositive integers), it will be converted to the sequence \code{1:v}. All atomic types are supported (See \code{?is.atomic}).}
  \item{numGroups}{An Integer. The number of groups that the vector will be partitioned into. Must divide the length of v (if v is a vector) or v (if v is a scalar).
}
}
\value{
A numerical value representing the total number of partitions of groups of equal size.
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
}
