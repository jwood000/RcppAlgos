\name{comboCount}
\alias{comboCount}
\alias{comboCount.default}
\alias{comboCount.table}
\alias{comboCount.list}
\alias{permuteCount}
\alias{permuteCount.default}
\alias{permuteCount.table}
\alias{permuteCount.list}
\title{
Number of combinations/permutations
}
\description{
Calculate the number of combinations/permutations of a vector chosen \eqn{m} at a time with or without replacement. Additionally, these functions can calculate the number of combinations/permutations of multisets.
}
\usage{
comboCount(v, m = NULL, ...)
permuteCount(v, m = NULL, ...)

\method{comboCount}{default}(v, m = NULL, repetition = FALSE, freqs = NULL, ...)
\method{permuteCount}{default}(v, m = NULL, repetition = FALSE, freqs = NULL, ...)

\method{comboCount}{table}(v, m = NULL, ...)
\method{permuteCount}{table}(v, m = NULL, ...)
}
\arguments{
  \item{v}{Source vector. If \code{v} is a positive integer, it will be converted to the sequence \code{1:v}. If \code{v} is a negative integer, it will be converted to the sequence \code{v:-1}. All atomic types are supported (See \code{\link{is.atomic}}).}
  \item{m}{Number of elements to choose.  If \code{repetition = TRUE} or \code{freqs} is utilized, \code{m} can exceed the length of \code{v}. If \code{m = NULL}, the length will default to \code{length(v)} or \code{sum(freqs)}.}
  \item{\dots}{Further arguments passed to methods.}
  \item{repetition}{Logical value indicating whether combinations/permutations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all combinations/permutations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
}
\value{
A numerical value representing the total number of combinations/permutations.
}
\note{
When the number of results exceeds \eqn{2^{53} - 1}{2^53 - 1}, a number of class \code{bigz} is returned.
}
\seealso{
\code{\link{comboGeneral}}, \code{\link{permuteGeneral}}
}
\examples{
## Same interface as the respective "general" functions:
## i.e. comboGeneral & permuteGeneral

permuteCount(-5)
permuteCount(5)
comboCount(25, 12)
permuteCount(15, 7, TRUE)
comboCount(25, 12, freqs = rep(2, 25))

## Return object of class 'bigz'
comboCount(250, 15, freqs = rep(2, 250))
}
\keyword{combinations}
\keyword{permutations}
\keyword{combinatorics}
