\name{comboCount}
\alias{comboCount}
\alias{permuteCount}
\title{
Number of combinations/permutations
}
\description{
Calculate the number of combinations/permutations of a vector chosen \eqn{m} at a time with or without replacement. Additionally, these functions can calculate the number of combinations/permutations of multisets.
}
\usage{
comboCount(v, m = NULL, repetition = FALSE, freqs = NULL)

permuteCount(v, m = NULL, repetition = FALSE, freqs = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is an integer, it will be converted to the sequence \code{1:v}.}
  \item{m}{Number of elements to choose. If \code{repetition = TRUE}, \code{m} can exceed the length of \code{v}. The default is \code{NULL}.}
  \item{repetition}{Logical value indicating whether combinations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all permutations/combinations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
}
\value{
A numerical value.
}
\note{
Accurate up to \eqn{2^{53} - 1}{2^53 - 1}.
}
\seealso{
\code{\link{comboGeneral}}, \code{\link{permuteGeneral}}
}
\examples{
## Same interface as the respective "general" functions:
## i.e. comboGeneral & permuteGeneral

comboCount(25, 12)

permuteCount(15, 7, TRUE)

comboCount(25, 12, freqs = rep(2, 25))
}
