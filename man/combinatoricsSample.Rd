\name{comboSample}
\alias{comboSample}
\alias{permuteSample}
\title{
Sample combinations/permutations
}
\description{
Generate a specific (lexicographically) or random sample of combinations/permutations.
}
\usage{
comboSample(v, m = NULL, repetition = FALSE, 
            freqs = NULL, n = NULL, sampleVec = NULL)

permuteSample(v, m = NULL, repetition = FALSE,
              freqs = NULL, n = NULL, sampleVec = NULL)
}
\arguments{
\item{v}{Source vector. If \code{v} is an integer, it will be converted to the sequence \code{1:v}.}
  \item{m}{Number of elements to choose. If \code{repetition = TRUE}, \code{m} can exceed the length of \code{v}. The default is \code{NULL}.}
  \item{repetition}{Logical value indicating whether combinations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all permutations/combinations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{n}{Number of combinations/permutations to return. The default is \code{NULL}.
}
  \item{sampleVec}{A vector of numbers representing the lexicographical combination/permutations to return.
}
}
\details{
These algorithms rely on efficiently generating the \emph{nth} lexicographical combination/permutation (sometimes called the \href{https://rosettacode.org/wiki/Permutations/Rank_of_a_permutation}{rank}).
}
\value{
In general, a matrix is returned with each row containing a vector of length \eqn{m}. If \code{m} isn't supplied and \code{freqs} is given, a matrix is returned with each row containing a vector of length \code{sum(freqs)}.

}
\references{
\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical order}
}
\author{
Joseph Wood
}
\note{
\itemize{
\item \code{n} and \code{sampleVec} cannot both be \code{NULL}.
\item Factor vectors are accepted. Class and level attributes are preserved.
}
}

\examples{
set.seed(11)
comboSample(30, 15, TRUE, n = 10)

set.seed(2)
permuteSample(15, 10, freqs = sample(1:2, 15, TRUE), 
              sampleVec = c(1, 10^2, 10^5, 10^8, 10^11))
              
all.equal(comboSample(10, 5, 
            sampleVec = 1:comboCount(10, 5)),
         comboGeneral(10, 5))
}